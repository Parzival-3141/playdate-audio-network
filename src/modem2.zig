const std = @import("std");

pub const OscillatorRates = extern struct {
    clock: [2]f32,
    header: [3]f32,
    payload: [8]f32,

    const osc_count = 13;

    /// Returns a set of oscillator phase rates (inverse frequencies) that are optimal for
    /// Goertzel detection with N DFT terms at the given sample rate.
    ///
    /// One symbol requires 13 unique frequencies, in increasing order:
    /// - Two for the clock oscillator, alternating low/high (FSK)
    /// - Three for the header oscillator, disconnected/connected/payload (FSK)
    /// - Eight for each payload bit oscillator, in least to most signifcant order (ASK)
    pub fn init(N: u16, sample_rate: f32) OscillatorRates {
        const base_freq = sample_rate / @as(f32, @floatFromInt(N));
        if (base_freq < 20) {
            @compileError("modem signal contains a frequency below 20 Hz; choose a smaller value for N");
        }
        if (base_freq * osc_count > 20_000) {
            @compileError("modem signal contains a frequency over 20 kHz; choose a larger value for N");
        }

        var rates: [osc_count]f32 = undefined;
        for (&rates, 0..) |*f, i| {
            f.* = (base_freq * (i + 1)) / sample_rate;
        }

        return @bitCast(rates);
    }
};

/// Creates a type that modulates data as an audio signal.
/// The signal uses a combination of frequency-shift keying (FSK) and amplitude-shifk keying (ASK).
/// Sine waves are generated using a lookup table of the given length.
/// 64 is a reasonable table length.
pub fn Modulator(
    comptime N: u16,
    comptime sample_rate: comptime_float,
    comptime baud: comptime_float,
    comptime sine_table_len: u32,
) type {
    const rates = OscillatorRates.init(N, sample_rate);

    return struct {
        const Mod = @This();

        const N_float: comptime_float = @floatFromInt(N);
        const window_len_float = sample_rate / baud;
        pub const window_len: u16 = @intFromFloat(window_len_float);
        comptime {
            if (window_len_float != @as(comptime_float, @floatFromInt(window_len))) {
                @compileError("baud must be integer divisible by sample rate");
            }
        }

        pub const sine_scaling = 0.2;
        pub const sine_table = create_sine_table(sine_table_len);

        pub const Symbol = union(enum) {
            disconnected,
            connected,
            payload: u8,
        };

        clock_osc: Oscillator,
        clock_osc_fm: Oscillator,
        header_osc: Oscillator,
        payload_oscs: [8]Oscillator,
        last_header: u2,
        last_payload: u8,

        pub fn init() Mod {
            return Mod{
                .clock_osc = .{},
                .clock_osc_fm = .{
                    // Initial phase is offset such that the clock freqs are strongest
                    // at the center of the the significant N samples of the window.
                    .phase = 1.0 - (N_float / 2.0) / window_len_float,
                },
                .header_osc = .{},
                .payload_oscs = .{.{}} ** 8,
                .last_header = 0,
                .last_payload = 0,
            };
        }

        /// Generate signal for the value (or no-op) into buf.
        /// Buf must be large enough to store one symbol period.
        pub fn modulate(mod: *Mod, symbol: Symbol, buf: []f32) void {
            for (buf[0..window_len]) |*out| {
                // Get clock FM value, which oscillates between the two clock frequencies.
                // TODO: This should be precalclated in OscillatorRates.
                const diff = rates.clock[1] - rates.clock[0];
                const center = rates.clock[0] + (diff / 2);
                var carrier_rate = mod.clock_osc_fm.generate(&sine_table, baud / sample_rate / 2.0); // [-1, 1]
                carrier_rate *= diff; // [-diff, diff]
                carrier_rate += center; // [center-diff, center+diff] = [low, high]

                // Clear buf on the first pass rather than summing existing signal,
                // which is expected to be undefined.
                out.* = mod.clock_osc.generate(&sine_table, carrier_rate) * sine_scaling;
            }

            // The number of samples used to smoothly ramp from an old frequency or amplitude to a new one.
            // This transition helps avoid creating bursts of energy in higher frequencies (noise).
            const transition_samples = window_len - N;

            {
                var start: u32 = 0;
                var len: u32 = window_len;

                const header = @intFromEnum(symbol);
                if (header != mod.last_header) {
                    // Transition from old frequency to new
                    for (buf[0..transition_samples], 0..) |*out, i| {
                        const frac = @as(f32, @floatFromInt(i)) / transition_samples;
                        const freq = cosine_interpolate(
                            &sine_table,
                            rates.header[mod.last_header],
                            rates.header[header],
                            frac,
                        );
                        var samp = mod.header_osc.generate(&sine_table, freq);
                        out.* += samp * sine_scaling;
                    }

                    start = transition_samples;
                    len = N;
                }
                // Hold frequency for remaining samples
                for (buf[start..][0..len]) |*out| {
                    var samp = mod.header_osc.generate(&sine_table, rates.header[header]);
                    out.* += samp * sine_scaling;
                }

                mod.last_header = header;
            }

            const byte = switch (symbol) {
                .payload => |b| b,
                else => 0,
            };

            for (0..8) |bit| {
                const mask: u8 = @as(u8, 1) << @intCast(bit);
                const bit_set = byte & mask != 0;
                const last_bit_set = mod.last_payload & mask != 0;

                // No output if holding a zero amplitude
                if (!bit_set and !last_bit_set) {
                    continue;
                }

                var start: u32 = 0;
                var len: u32 = window_len;

                if (bit_set != last_bit_set) {
                    // Transition from old amplitude to new
                    for (buf[0..transition_samples], 0..) |*out, i| {
                        const frac = @as(f32, @floatFromInt(i)) / transition_samples;
                        const amp = cosine_interpolate(
                            &sine_table,
                            @floatFromInt(@intFromBool(last_bit_set)),
                            @floatFromInt(@intFromBool(bit_set)),
                            frac,
                        );
                        var samp = mod.payload_oscs[bit].generate(&sine_table, rates.payload[bit]);
                        out.* += samp * amp * (sine_scaling / 2.0);
                    }

                    start = transition_samples;
                    len = N;
                }

                // Hold amplitude for remaining samples
                if (!bit_set) {
                    continue;
                }
                for (buf[start..][0..len]) |*out| {
                    var samp = mod.payload_oscs[bit].generate(&sine_table, rates.payload[bit]);
                    out.* += samp * (sine_scaling / 2.0);
                }
            }

            mod.last_payload = byte;
        }
    };
}

test Modulator {
    const M = Modulator(31, 44100, 882, 64);
    var m = M.init();
    var buf: [100]f32 = undefined;
    m.modulate(.disconnected, &buf);
    m.modulate(.disconnected, &buf);
    m.modulate(.connected, &buf);
    m.modulate(.{ .payload = 'h' }, &buf);
    m.modulate(.{ .payload = 'i' }, &buf);
    m.modulate(.connected, &buf);
}

/// Creates a sinusoidal lookup table of the given length.
/// Holds amplitudes in range [-1, 1] for one cycle of a sine.
/// A length of 64 is fine for many practical purposes.
/// Lookup table is scaled ahead of time to compensate for expected amplitude summation at runtime.
pub fn create_sine_table(comptime len: u32) [len]f32 {
    if (len == 0) {
        @compileError("sine lookup table length must be > 0");
    }

    var vals: [len]f32 = undefined;
    for (&vals, 0..) |*v, i| {
        var frac: f64 = @floatFromInt(i);
        frac /= vals.len;
        v.* = @floatCast(@sin(frac * std.math.pi * 2));
    }

    return vals;
}

/// Returns sin(phase * 2pi) using the lookup table with linear interpolation.
/// Phase must be in range [0, 1).
pub fn sine_value(sine_table: []const f32, phase: f32) f32 {
    std.debug.assert(phase >= 0 and phase < 1);

    const index_float = phase * @as(f32, @floatFromInt(sine_table.len));
    const index: u32 = @intFromFloat(index_float);
    const frac = index_float - @as(f32, @floatFromInt(index));

    const y1 = sine_table[index % sine_table.len];
    const y2 = sine_table[(index + 1) % sine_table.len];

    return std.math.lerp(y1, y2, frac);
}

/// Return a value that is a fraction of the way between `from` and `to`, where the transition is
/// curved by a half cosine rather than linear. Asserts frac is in range [0, 1).
pub fn cosine_interpolate(sine_table: []const f32, from: f32, to: f32, frac: f32) f32 {
    std.debug.assert(frac >= 0 and frac < 1);

    var p = frac / 2;
    p += 0.25;

    var y = sine_value(sine_table, p); // [1, -1]
    y /= 2; // [0.5, -0.5]
    y += 0.5; // [1, 0]
    y = 1 - y; // [0, 1]

    const diff = to - from;
    return from + (y * diff);
}

test cosine_interpolate {
    const sine_table = create_sine_table(128);
    for (0..100) |i| {
        const frac = @as(f32, @floatFromInt(i)) / 100.0;
        const expected = @cos((1.0 + frac) * std.math.pi) / 2.0 + 0.5;
        const res = cosine_interpolate(&sine_table, 0, 1, frac);
        try std.testing.expectApproxEqAbs(expected, res, 0.0005);
    }
}

pub fn Demodulator(comptime N: u16, comptime sample_rate: f32) type {
    const freqs = SymbolFrequencies.init(N, sample_rate);

    return struct {
        pub const Demod = @This();

        pub const Error = error{NoSignal};

        pub const State = enum {
            start,
            clock_low,
            clock_high,
        };

        // last_clock: bool,
        state: State,
        // buf: [N]f32,
        // buf_index: u16,

        pub fn init() Demod {
            return .{
                // .last_clock = false,
                .state = .start,
                // .buf = &[0]f32{},
                // .buf_index = 0,
            };
        }

        /// Get the ratio of the higher value to the lower value.
        /// Bounds are set to prevent NaN and Inf.
        /// Also returns the higher of the two values.
        fn magnitude_ratio(a: f32, b: f32) struct { f32, u1 } {
            std.debug.assert(a >= 0);
            std.debug.assert(b >= 0);

            var lower = a;
            var higher = b;
            var sel: u1 = 1;
            if (lower > higher) {
                lower = b;
                higher = a;
                sel = 0;
            }

            lower = std.math.clamp(lower, 0.00001, 100.0);
            higher = std.math.clamp(higher, 0.00001, 100.0);

            return .{ higher / lower, sel };
        }

        pub fn demodulate(demod: *Demod, signal: []const f32) Error!?u8 {
            std.log.debug("----", .{});

            var high_mag: f32 = undefined;
            {
                const clk0 = goertzel(N, signal, freqs.clock[0]);
                const clk1 = goertzel(N, signal, freqs.clock[1]);
                // TODO: upgrade zig and destructure
                const cmp = magnitude_ratio(clk0, clk1);
                const ratio = cmp[0];
                const higher = cmp[1];

                std.log.debug("clk: {d} | {d} | {d}", .{ clk0, clk1, ratio });

                // TODO: resync here
                if (ratio < 100) {
                    demod.state = .start;
                    return error.NoSignal;
                }
                const high = higher == 1;
                high_mag = ([2]f32{ clk0, clk1 })[higher];

                // demod.last_clock = high;

                // Same clock value twice in a row results in a re-sync.
                switch (demod.state) {
                    .start => {
                        demod.state = if (high) .clock_high else .clock_low;
                    },
                    .clock_low => {
                        if (high) {
                            demod.state = .clock_high;
                        } else {
                            demod.state = .start;
                            return error.NoSignal;
                        }
                    },
                    .clock_high => {
                        if (high) {
                            demod.state = .start;
                            return error.NoSignal;
                        } else {
                            demod.state = .clock_low;
                        }
                    },
                }
            }

            const data = blk: {
                const data0 = goertzel(N, signal, freqs.data[0]);
                const data1 = goertzel(N, signal, freqs.data[1]);
                // TODO: upgrade zig and destructure
                const cmp = magnitude_ratio(data0, data1);
                const ratio = cmp[0];
                const higher = cmp[1];

                std.log.debug("data: {d} | {d} | {d}", .{ data0, data1, ratio });

                // TODO: resync here
                if (ratio < 100) {
                    demod.state = .start;
                    return error.NoSignal;
                }
                const high = higher == 1;

                break :blk high;
            };

            if (!data) {
                return null;
            }

            var byte: u8 = 0;
            for (freqs.payload, 0..) |freq, i| {
                const mag = goertzel(N, signal, freq);
                const target_mag = high_mag / 20;
                std.log.debug("freq {d} | mag {d} | target {d}", .{ i, mag, target_mag });
                const val: u8 = @intFromBool(mag > target_mag); // TODO
                const mask: u8 = val << @intCast(i);
                byte |= mask;
            }

            return byte;
        }
    };
}

/// Calculates the power of the given frequency within the input sample slice.
/// N is the number of DFT terms. Buffer must have at least N samples.
/// normalized_frequency is target frequency divided by sample rate.
/// Frequency detection is strongest when frequency is an integer multiple of (1/N * sample_rate).
///
/// Reference: https://en.wikipedia.org/wiki/Goertzel_algorithm#Power-spectrum_terms
pub fn goertzel(N: u16, buf: []const f32, normalized_frequency: f32) f32 {
    const coeff = 2 * @cos(2.0 * std.math.pi * normalized_frequency);

    var s_prev: f32 = 0;
    var s_prev2: f32 = 0;

    for (0..N) |i| {
        const s = buf[i] + (coeff * s_prev) - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }

    // TODO: scale power based on N
    return (s_prev2 * s_prev2) + (s_prev * s_prev) - (coeff * s_prev * s_prev2);
}

/// Sine wave oscillator that uses a lookup table with linear interpolation.
pub const Oscillator = struct {
    /// Current phase of the sine wave in range [0, 1).
    phase: f32 = 0,

    var called: usize = 0;

    /// Produce the next oscillator amplitude and advance phase by phase_incr.
    /// Phase increment can be calculated with (freq / sample_rate);
    /// it is provided as a parameter instead of frequency because for this package because
    /// it is more effecient to precalculate the phase_incrs of a finite set of frequencies.
    ///
    /// phase_incr is asserted to be positive.
    pub fn generate(osc: *Oscillator, comptime sine_table: []const f32, phase_incr: f32) f32 {
        std.debug.assert(phase_incr >= 0);
        std.debug.assert(phase_incr < 0.5);
        std.debug.assert(osc.phase <= 1 and osc.phase >= 0);

        const index_float = osc.phase * sine_table.len;
        const index: u32 = @intFromFloat(index_float);
        const frac = index_float - @as(f32, @floatFromInt(index));

        const y1 = sine_table[index % sine_table.len];
        const y2 = sine_table[(index + 1) % sine_table.len];
        const out = y1 * (1 - frac) + (y2 * frac);

        osc.phase += phase_incr;
        while (osc.phase >= 1) osc.phase -= 1;

        return out;
    }
};

test "modulate and demodulate symbols" {
    const N = 30;

    var mod = Modulator(N, 44_100, 64).init();
    var demod = Demodulator(N, 44_100).init();

    const payloads = [_]?u8{
        null, null, null, null,
        'z',  'i',  'g',  '.',
        '.',  '.',  ' ',  null,
        null, null, null, null,
        'H',  'I',  ',',  ' ',
        'M',  'O',  'M',  '!',
        null, null, null, null,
    };

    var signal = [_]f32{0} ** (N * payloads.len);

    for (payloads, 0..) |p, i| {
        mod.modulate(p, signal[i * N ..]);
    }

    for (payloads, 0..) |p, i| {
        const sym = try demod.demodulate(signal[N * i ..]);

        // TODO: Use testing log instead of panicking to see which iteration failed.
        // https://github.com/ziglang/zig/issues/5738
        std.testing.expectEqual(p == null, sym == null) catch {
            std.debug.panic("payload null mismatch on iter {}", .{i});
        };
        std.testing.expectEqual(p, sym) catch {
            std.debug.panic("symbol mismatch on iter {}: expected {?c} ({?x}), got {?c} ({?x})", .{ i, p, p, sym, sym });
        };
    }
}
