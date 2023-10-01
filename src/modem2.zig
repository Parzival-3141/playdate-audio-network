const std = @import("std");

pub const SymbolFrequencies = extern struct {
    clock: [2]f32,
    data: [2]f32,
    payload: [8]f32,

    /// Returns a set of frequencies that are optimal for Goertzel detection
    /// with N DFT terms at the given sample rate.
    /// Each frequency is normalized as a fraction of sample rate, which is also the
    /// oscillator's per-rample phase increment.
    ///
    /// One symbol requires 12 unique frequencies, in increasing order:
    /// - Two for the clock bit oscillator, alternating low/high (FSK)
    /// - Two for the data bit oscillator, no/yes (FSK)
    /// - Eight for each payload bit oscillator, in least to most signifcant order (ASK)
    pub fn init(N: u16, sample_rate: f32) SymbolFrequencies {
        const base_freq = sample_rate / @as(f32, @floatFromInt(N));
        if (base_freq < 20) {
            @compileError("modem signal contains a frequency below 20 Hz; choose a smaller value for N");
        }
        if (base_freq * 12 > 20_000) {
            @compileError("modem signal contains a frequency over 20 kHz; choose a larger value for N");
        }

        var freqs: [12]f32 = undefined;
        for (&freqs, 0..) |*f, i| {
            f.* = (base_freq * (i + 1)) / sample_rate;
        }

        return @bitCast(freqs);
    }
};

/// Creates a type that modulates data as an audio signal.
/// The signal uses a combination of frequency-shift keying (FSK) and amplitude-shifk keying (ASK).
/// Sine waves are generated using a lookup table of the given length.
/// 64 is a reasonable table length.
pub fn Modulator(comptime N: u16, comptime sample_rate: f32, comptime sine_table_len: u32) type {
    const freqs = SymbolFrequencies.init(N, sample_rate);

    return struct {
        const Mod = @This();

        pub const sine_table = create_sine_table(sine_table_len, 0.25);

        clock: bool,
        clock_osc: Oscillator,
        data_osc: Oscillator,
        payload_oscs: [8]Oscillator,

        pub fn init() Mod {
            return Mod{
                .clock = false,
                .clock_osc = .{},
                .data_osc = .{},
                .payload_oscs = .{.{}} ** 8,
            };
        }

        /// Generate signal for the value (or no-op) into buf.
        /// Buf must be large enough to store one symbol period.
        pub fn modulate(mod: *Mod, value: ?u8, buf: []f32) void {
            // TODO: Check performance of summing one sample at a time vs.
            // adding all clock bit samples, then data bit samples, then payload samples.
            // The latter is probably more optimized.

            for (buf[0..N]) |*out| {
                // Clear buf on the first pass rather than summing existing signal,
                // which is expected to be undefined.
                out.* = mod.clock_osc.generate(&sine_table, freqs.clock[@intFromBool(mod.clock)]);
            }

            const data_bit = value != null;
            for (buf[0..N]) |*out| {
                out.* += mod.data_osc.generate(&sine_table, freqs.data[@intFromBool(data_bit)]);
            }

            if (value) |byte| {
                for (0..8) |bit| {
                    const mask: u8 = @as(u8, 1) << @intCast(bit);
                    if (byte & mask == 0) {
                        continue;
                    }

                    for (buf[0..N], 0..N) |*out, i| {
                        _ = i;
                        var samp = mod.payload_oscs[bit].generate(&sine_table, freqs.payload[bit]);

                        // // Apply triangular window to output signal so that sudden oscillator
                        // // starts and stops do not produce sidebands in the DAC and ADC.
                        // //
                        // // TODO: Explore other windowing functions that preserve more energy
                        // // of the desired frequency. Maybe Welch window?

                        // const i2_float: f32 = @floatFromInt(i * 2);
                        // const N_float: f32 = @floatFromInt(N);
                        // const frac = i2_float / (N_float); // Progress through the window, [0, 2]
                        // if (frac <= 1) {
                        //     samp *= frac;
                        // } else {
                        //     samp *= (2 - frac);
                        // }

                        out.* += samp * 0.5;
                    }
                }
            }

            mod.clock = !mod.clock;
        }
    };
}

/// Creates a sinusoidal lookup table of the given length.
/// Holds amplitudes in range [-1, 1] for one cycle of a sine.
/// A length of 64 is fine for many practical purposes.
pub fn create_sine_table(comptime len: u32, scale: f32) [len]f32 {
    if (len == 0) {
        @compileError("sine lookup table length must be > 0");
    }
    if (scale < 0 or scale > 1) {
        std.debug.panic("scale must be between 0 and 1, got {}", .{scale});
    }

    var vals: [len]f32 = undefined;
    for (&vals, 0..) |*v, i| {
        v.* = @floatCast(@sin(i / @as(f64, vals.len) * std.math.pi * 2));
        v.* *= scale;
    }

    return vals;
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

        pub fn demodulate(demod: *Demod, signal: []const f32) Error!?u8 {
            {
                const clk0 = goertzel(N, signal, freqs.clock[0]);
                const clk1 = goertzel(N, signal, freqs.clock[1]);
                const high = clk1 > clk0;
                const diff = if (high) clk1 - clk0 else clk0 - clk1;
                // TODO: if diff is below threshold, re-sync. if diff still low, error.NoSignal
                // TODO: what is a reasonable diff?
                if (diff < 0.5) {
                    demod.state = .start;
                    return error.NoSignal;
                }

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
                const high = data1 > data0;
                const diff = if (high) data1 - data0 else data0 - data1;

                // TODO: what is a reasonable diff?
                if (diff < 0.5) {
                    demod.state = .start;
                    return error.NoSignal;
                }

                break :blk high;
            };

            if (!data) {
                return null;
            }

            var byte: u8 = 0;
            for (freqs.payload, 0..) |freq, i| {
                const mag = goertzel(N, signal, freq);
                const val: u8 = @intFromBool(mag > 0.5); // TODO
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
        std.testing.expectEqual(p == null, sym == null) catch {
            std.debug.panic("payload null mismatch on iter {}", .{i});
        };
        std.testing.expectEqual(p, sym) catch {
            std.debug.panic("symbol mismatch on iter {}: expected {?c} ({?x}), got {?c} ({?x})", .{ i, p, p, sym, sym });
        };
    }
}
