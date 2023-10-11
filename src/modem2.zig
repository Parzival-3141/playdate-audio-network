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
    /// - Three for the header oscillator, waiting/ready/payload (FSK)
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

pub const Symbol = union(enum) {
    waiting,
    ready,
    payload: u8,
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
        const symbol_len_float = sample_rate / baud;
        pub const symbol_len: u16 = @intFromFloat(symbol_len_float);
        comptime {
            if (symbol_len_float != @as(comptime_float, @floatFromInt(symbol_len))) {
                @compileError("baud must be integer divisible by sample rate");
            }
        }

        pub const sine_scaling = 0.2;
        pub const sine_table = create_sine_table(sine_table_len);

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
                    .phase = 1.0 - (N_float / 2.0) / symbol_len_float,
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
            for (buf[0..symbol_len]) |*out| {
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
            const transition_samples = symbol_len - N;

            {
                var start: u32 = 0;
                var len: u32 = symbol_len;

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
                var len: u32 = symbol_len;

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
    m.modulate(.waiting, &buf);
    m.modulate(.waiting, &buf);
    m.modulate(.ready, &buf);
    m.modulate(.{ .payload = 'h' }, &buf);
    m.modulate(.{ .payload = 'i' }, &buf);
    m.modulate(.ready, &buf);
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

pub fn Demodulator(
    comptime N: u16,
    comptime sample_rate: comptime_float,
    comptime baud: comptime_float,
) type {
    if (baud <= 0) {
        @compileError("baud must be greater than 0");
    }

    const rates = OscillatorRates.init(N, sample_rate);

    return struct {
        const Demod = @This();

        const N_float: comptime_float = @floatFromInt(N);

        const symbol_len_float = sample_rate / baud;
        pub const symbol_len: u16 = @intFromFloat(symbol_len_float);
        comptime {
            if (symbol_len_float != @as(comptime_float, @floatFromInt(symbol_len))) {
                @compileError("baud must be integer divisible by sample rate");
            }
        }

        const padding_count = symbol_len - N;
        const lookaround_count = 2;

        comptime {
            std.debug.assert(lookaround_count < N);
            std.debug.assert(lookaround_count < padding_count);
        }
        const near_sync_samples_needed = N + (2 * lookaround_count);
        const full_sync_samples_needed = symbol_len + N;

        pub const ClockState = union(enum) {
            start,
            set: u1,
        };

        clock_state: ClockState,
        near_sync_countdown: u16,
        /// Buffer of incoming samples taken from the incoming signal stream.
        /// Only used when the slice passed to demodulate() does not have enough samples to support
        /// the given analysis (normally N, but more for full or partial sync).
        buffer: [full_sync_samples_needed]f32,
        /// How many values are currently in the buffer.
        buffer_items: u16,
        /// How many samples should be dropped before acquiring the next meaningful slice.
        buffer_skip_samples: u16,

        pub fn init() Demod {
            return .{
                .clock_state = .start,
                .near_sync_countdown = 10,
                .buffer = undefined,
                .buffer_items = 0,
                .buffer_skip_samples = 0,
            };
        }

        pub const Result = struct {
            /// Number of samples read.
            /// Caller should advance signal cursor by this amount after calling.
            usize,
            /// Status of analysis.
            union(enum) {
                /// The incoming signal did not yield a strong signal,
                /// or the detected clock was an unexpected value.
                /// On this status, state is reset and the next call to demodulate()
                /// will invoke a full re-sync.
                disconnected,
                /// Signal was in tact.
                symbol: Symbol,
            },
        };

        /// Reads as many samples as needed from signal to decode a single symbol,
        /// buffering internally for future calls if needed. Once enough samples are read,
        /// returns either the decoded symbol or .disconnected if the remote signal is not
        /// strong or clear enough.
        ///
        /// Returns null if more samples are needed to decode the next symbol, in which case
        /// the end of `signal` was reached.
        pub fn demodulate(d: *Demod, signal: []const f32) ?Result {
            std.log.debug("demodulating symbol", .{});

            const clock_res = switch (d.clock_state) {
                .start => res: {
                    break :res d.check_clock_full_sync(signal, 6) orelse return null;
                },
                .set => res: {
                    if (d.near_sync_countdown == 0) {
                        d.near_sync_countdown = 10;
                        break :res d.check_clock_near_sync(signal) orelse return null;
                    } else {
                        d.near_sync_countdown -= 1;
                        break :res d.check_clock_no_sync(signal) orelse return null;
                    }
                },
            };

            const read, const status = clock_res;
            const analysis_slice = switch (status) {
                .weak_signal, .wrong_clock => return .{ read, .disconnected },
                .ok => |slice| slice,
            };

            {
                const hdr0 = goertzel(N, analysis_slice, rates.header[0]);
                const hdr1 = goertzel(N, analysis_slice, rates.header[1]);
                const hdr2 = goertzel(N, analysis_slice, rates.header[2]);
                const sel, const ratio = select_magnitude(&.{ hdr0, hdr1, hdr2 });
                if (ratio < 1000) {
                    d.* = init();
                    return .{ read, .disconnected };
                }

                std.log.debug("hdr:  {d} | {d} | {d}", .{ hdr0, hdr1, hdr2 });

                switch (@as(u2, @intCast(sel))) {
                    0 => return .{ read, .{ .symbol = .waiting } },
                    1 => return .{ read, .{ .symbol = .ready } },
                    2 => {}, // below
                    3 => unreachable,
                }
            }

            var byte: u8 = 0;
            for (rates.payload, 0..) |rate, bit| {
                const mag = goertzel(N, analysis_slice, rate);
                const target_mag = 2.0; // TODO
                const val: u8 = @intFromBool(mag > target_mag); // TODO
                const mask: u8 = val << @intCast(bit);
                byte |= mask;

                std.log.debug("bit  {d} | mag {d} | target {d}", .{ bit, mag, target_mag });
            }

            return .{ read, .{ .symbol = .{ .payload = byte } } };
        }

        /// When successful, returns the amount read and a slice of `needed` samples.
        /// If signal is long enough and internal buffering is not already in progress,
        /// the signal slice is used directly to avoid copying.
        /// Otherwise, buffers incoming signal until `needed` samples has been collected.
        /// Until enough samples are ready, returns null.
        fn collect_signal(d: *Demod, signal: []const f32, needed: usize) ?struct { usize, []const f32 } {
            if (d.buffer_skip_samples > signal.len) {
                d.buffer_skip_samples -= @intCast(signal.len);
                return null;
            }

            const sig = signal[d.buffer_skip_samples..];
            d.buffer_skip_samples = 0;

            // Happy path: no buffering was needed.
            if (d.buffer_items == 0 and sig.len >= needed) {
                return .{ needed, sig[0..needed] };
            }

            // Otherwise, either the signal slice is too short or we need to append to
            // an in-progress buffering.
            const copy_count: u16 = @intCast(needed - d.buffer_items);
            if (copy_count < sig.len) {
                @memcpy(d.buffer[d.buffer_items..][0..copy_count], sig[0..copy_count]);
                d.buffer_items += copy_count;
                return .{ copy_count, d.buffer[0..needed] };
            } else {
                @memcpy(d.buffer[d.buffer_items..][0..sig.len], sig);
                d.buffer_items += @intCast(sig.len);
                return null;
            }
        }

        test collect_signal {
            var d = Demod.init();
            const static_signal: [full_sync_samples_needed]f32 = undefined;

            const Step = struct {
                skip: ?u16 = null,
                sig: []const f32,
                needed: u16,
                res: ?struct { usize, []const f32 },
                buffer_items: u16,
                reset: bool = false,
            };
            const steps: []const Step = &.{
                .{
                    .sig = &static_signal,
                    .needed = 0,
                    .res = .{ 0, static_signal[0..0] },
                    .buffer_items = 0,
                },
                .{
                    .sig = &static_signal,
                    .needed = 1,
                    .res = .{ 1, static_signal[0..1] },
                    .buffer_items = 0,
                },
                .{
                    .sig = static_signal[0..symbol_len],
                    .needed = symbol_len,
                    .res = .{ symbol_len, static_signal[0..symbol_len] },
                    .buffer_items = 0,
                },
                .{
                    .sig = static_signal[0..symbol_len],
                    .needed = symbol_len + 1,
                    .res = null,
                    .buffer_items = symbol_len,
                },
                .{
                    .sig = static_signal[0..symbol_len],
                    .needed = symbol_len + 1,
                    .res = .{ 1, d.buffer[0 .. symbol_len + 1] },
                    .buffer_items = symbol_len + 1,
                },
                .{ // part 1
                    .reset = true,
                    .skip = full_sync_samples_needed / 2,
                    .sig = &static_signal,
                    .needed = full_sync_samples_needed,
                    .res = null,
                    .buffer_items = static_signal.len - (full_sync_samples_needed / 2),
                },
                .{ // part 2
                    .sig = &static_signal,
                    .needed = full_sync_samples_needed,
                    .res = .{ full_sync_samples_needed / 2, d.buffer[0..full_sync_samples_needed] },
                    .buffer_items = full_sync_samples_needed,
                },
            };

            for (steps) |step| {
                if (step.reset) {
                    d.buffer_items = 0;
                }
                if (step.skip) |amount| {
                    d.buffer_skip_samples = amount;
                }
                const res = d.collect_signal(step.sig, step.needed);
                try std.testing.expectEqual(step.res, res);
                try std.testing.expectEqual(step.buffer_items, d.buffer_items);
            }
        }

        const SyncResult = struct {
            usize,
            union(enum) {
                weak_signal,
                wrong_clock,
                ok: []const f32,
            },
        };

        /// Analyzes signal to try to find clock peak, corresponding to
        /// optimal analysis offset. Does a binary search for the highest magnitude
        /// ratio of one clock value to the other over a full symbol period.
        /// If the signal slice does not have enough samples for that check,
        /// buffers as many as possible for the next call.
        ///
        /// Whichever clock value is detected as strongest is set as the clock state,
        /// and successive windows are expected to alternate clock values.
        fn check_clock_full_sync(d: *Demod, signal: []const f32, max_iters: u8) ?SyncResult {
            const read, const sig = d.collect_signal(signal, full_sync_samples_needed) orelse return null;

            var center: f32 = symbol_len / 2;
            var distance: f32 = symbol_len / 4;
            var iter: u16 = 0;

            var clock_sel: u1 = undefined;
            var ratio: f32 = 1;
            var prev_clock_sel: u1 = undefined;
            var prev_ratio: f32 = 1;

            while (distance >= 0.5 and iter < max_iters) : ({
                iter += 1;
                distance /= 2;
                prev_clock_sel = clock_sel;
                prev_ratio = ratio;
            }) {
                const positions = [2]f32{
                    center - distance, // left
                    center + distance, // right
                };

                const left_chunk = sig[@intFromFloat(positions[0])..][0..N];
                const right_chunk = sig[@intFromFloat(positions[1])..][0..N];

                const mags = [4]f32{
                    goertzel(N, left_chunk, rates.clock[0]),
                    goertzel(N, left_chunk, rates.clock[1]),
                    goertzel(N, right_chunk, rates.clock[0]),
                    goertzel(N, right_chunk, rates.clock[1]),
                };

                const new_sel, const new_ratio = select_magnitude(&mags);
                if (new_ratio > prev_ratio) {
                    clock_sel = @as(u1, @intCast(new_sel % 2));
                    ratio = new_ratio;
                    center = positions[new_sel / 2];
                }
            }

            if (ratio < 1000) {
                d.* = Demod.init();
                return .{ read, .weak_signal };
            }
            d.clock_state = .{ .set = clock_sel };

            // Resolved analysis window.
            const start: usize = @intFromFloat(center);
            const slice = sig[start..][0..N];

            // TODO: It should be possible to "unread" by as much as lookaround_count
            // to avoid copying into the buffer so successive calls can always avoid
            // relying on the buffer when incoming signal is long enough.
            // As is, the buffer is always used.

            // The position of the next meaningful symbol portion, minus look-behind.
            const next_needed_sample_pos = start + symbol_len - lookaround_count;
            if (next_needed_sample_pos < d.buffer_items) {
                // Some of what we need for the next analysis (look-behind, N, look-ahead)
                // is already in the current slice. Move as many as possible to the start
                // of the buffer for the next run.
                const to_move = sig[next_needed_sample_pos..];
                for (to_move, d.buffer[0..to_move.len]) |src, *dst| {
                    dst.* = src;
                }
                d.buffer_skip_samples = 0;
                d.buffer_items = @intCast(to_move.len);
            } else {
                d.buffer_skip_samples = @intCast(next_needed_sample_pos - d.buffer_items);
                d.buffer_items = 0;
            }

            return .{ read, .{ .ok = slice } };
        }

        /// Analyzes signal to find clock peak and returns a slice of samples ready for symbol analysis.
        /// This should only be called some time after a successful full sync was performed.
        /// Differences with full sync:
        ///
        /// - Only searches a small, contiguous neighborhood around the current sync point.
        ///   Normal drift is expected to be slight, so this window is only a few samples on
        ///   either side of the current positioning.
        /// - Expects the clock value to oscillate correctly. Full sync just picks whichever
        ///   is strongest as the starting value.
        fn check_clock_near_sync(d: *Demod, signal: []const f32) ?SyncResult {
            const read, const sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;

            const new_clock = d.clock_state.set +% 1;
            const rate = rates.clock[new_clock];

            // Check a sliding window of N samples to the left and right of current offset.
            var start: u16 = 0;
            var best_mag: f32 = -1;
            for (0..(2 * lookaround_count)) |offset| {
                const slice = sig[offset..][0..N];
                const mag = goertzel(N, slice, rate);
                if (mag > best_mag) {
                    best_mag = mag;
                    start = @intCast(offset);
                } else break;
            }
            const analysis_slice = sig[start..][0..N];

            // Compare to opposite clock frequency strength to ensure the remote clock
            // is still oscillating correctly.
            const opposite_clock_mag = goertzel(N, analysis_slice, rates.clock[d.clock_state.set]);

            const sel, const ratio = select_magnitude(&.{ best_mag, opposite_clock_mag });
            if (ratio < 1000) {
                d.* = Demod.init();
                return .{ read, .weak_signal };
            }
            if (sel != 0) {
                d.* = Demod.init();
                return .{ read, .wrong_clock };
            }
            d.clock_state = .{ .set = new_clock };

            // TODO: It should be possible to "unread" by as much as lookaround_count
            // to avoid copying into the buffer so successive calls can always avoid
            // relying on the buffer when incoming signal is long enough.
            // As is, the buffer is always used.

            // The position of the next meaningful symbol portion, minus look-behind.
            const next_needed_sample_pos = start + symbol_len - lookaround_count;
            if (next_needed_sample_pos < d.buffer_items) {
                // Some of what we need for the next analysis (look-behind, N, look-ahead)
                // is already in the current slice. Move as many as possible to the start
                // of the buffer for the next run.
                const to_move = sig[next_needed_sample_pos..];
                for (to_move, d.buffer[0..to_move.len]) |src, *dst| {
                    dst.* = src;
                }
                d.buffer_skip_samples = 0;
                d.buffer_items = @intCast(to_move.len);
            } else {
                d.buffer_skip_samples = next_needed_sample_pos - d.buffer_items;
                d.buffer_items = 0;
            }

            return .{ read, .{ .ok = sig[start..][0..N] } };
        }

        /// Validates the clock of the incoming signal and returns a slice of samples ready for symbol analysis.
        fn check_clock_no_sync(d: *Demod, signal: []const f32) ?SyncResult {
            const read, var sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;
            sig = sig[lookaround_count..][0..N];

            const clk0 = goertzel(N, sig, rates.clock[0]);
            const clk1 = goertzel(N, sig, rates.clock[1]);
            const sel, const ratio = select_magnitude(&.{ clk0, clk1 });
            if (ratio < 1000) {
                d.* = Demod.init();
                return .{ read, .weak_signal };
            }
            if (sel == d.clock_state.set) {
                d.* = Demod.init();
                return .{ read, .wrong_clock };
            }
            d.clock_state = .{ .set = @intCast(sel) };

            const next_needed_sample_pos = symbol_len - lookaround_count;
            // TODO: ensure this is always true by requiring (2*lookaround) < padding
            std.debug.assert(next_needed_sample_pos >= d.buffer_items);
            d.buffer_skip_samples = next_needed_sample_pos - d.buffer_items;
            d.buffer_items = 0;

            return .{ read, .{ .ok = sig } };
        }

        test check_clock_no_sync {
            var d = Demod.init();
            d.clock_state = .{ .set = 1 };
            var sig: [near_sync_samples_needed]f32 = undefined;

            // Fill with strong wave at 1x baud (low clock frequency) and weak wave at 2x baud
            {
                var phase: f32 = 0;
                for (&sig) |*samp| {
                    samp.* = @sin(phase);
                    samp.* += (@sin(phase * 2) / 1200);
                    phase += rates.clock[0] * 2 * std.math.pi;
                }

                const res = d.check_clock_no_sync(&sig).?;
                const read, const status = res;
                try std.testing.expectEqual(near_sync_samples_needed, @intCast(read));
                try std.testing.expectEqual(@as([]const f32, sig[lookaround_count..][0..N]), status.ok);
                try std.testing.expectEqual(@as(u1, 0), d.clock_state.set);

                d.buffer_items = 0;
                d.buffer_skip_samples = 0;
            }

            // Expected clock has been flipped;
            // Fill with strong wave at 2x baud (low clock frequency) and weak wave at 1x baud
            {
                var phase: f32 = 0;
                for (&sig) |*samp| {
                    samp.* = @sin(phase);
                    samp.* += (@sin(phase / 2) / 1200);
                    phase += rates.clock[1] * 2 * std.math.pi;
                }

                const res = d.check_clock_no_sync(&sig).?;
                const read, const status = res;
                try std.testing.expectEqual(near_sync_samples_needed, @intCast(read));
                try std.testing.expectEqual(@as([]const f32, sig[lookaround_count..][0..N]), status.ok);
                try std.testing.expectEqual(@as(u1, 1), d.clock_state.set);

                d.buffer_items = 0;
                d.buffer_skip_samples = 0;
            }

            // Expected clock flipped again, so this should fail when given the previous signal
            {
                const res = d.check_clock_no_sync(&sig).?;
                const read, const status = res;
                try std.testing.expectEqual(near_sync_samples_needed, @intCast(read));
                try std.testing.expect(status == .wrong_clock);
            }
        }
    };
}

/// Returns the index of the highest value and the ratio of it to the second highest value.
/// Compared values are clamped to a reasonable range to prevent NaN and Inf.
fn select_magnitude(mags: []const f32) struct { u16, f32 } {
    std.debug.assert(mags.len > 0);

    var sel: u16 = 0;
    var max_mag: f32 = 0;
    var next_mag: f32 = 0;

    for (mags, 0..) |mag, i| {
        std.debug.assert(mag >= 0);
        if (mag > max_mag) {
            sel = @intCast(i);
            next_mag = max_mag;
            max_mag = mag;
        } else if (mag > next_mag) {
            next_mag = mag;
        }
    }

    next_mag = std.math.clamp(next_mag, 0.00001, 1000.0);
    max_mag = std.math.clamp(max_mag, 0.00001, 1000.0);

    return .{ sel, max_mag / next_mag };
}

test select_magnitude {
    const Case = struct {
        mags: []const f32,
        sel: u16,
        ratio: f32,
    };
    const cases: []const Case = &.{
        .{
            .mags = &.{ 10, 20 },
            .sel = 1,
            .ratio = 2.0,
        },
        .{
            .mags = &.{ 90, 20 },
            .sel = 0,
            .ratio = 4.5,
        },
        .{
            .mags = &.{ 10, 25, 0.5, 40, 20 },
            .sel = 3,
            .ratio = 1.6,
        },
    };
    for (cases) |c| {
        const sel, const ratio = select_magnitude(c.mags);
        try std.testing.expectEqual(c.sel, sel);
        try std.testing.expectApproxEqAbs(c.ratio, ratio, 0.001);
    }
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

    /// Produce the next oscillator amplitude and advance phase by phase_incr.
    /// Phase increment can be calculated with (freq / sample_rate);
    /// it is provided as a parameter instead of frequency because for this package because
    /// it is more effecient to precalculate the phase_incrs of a finite set of frequencies.
    ///
    /// phase_incr is asserted to be in the range [0, 0.5) to avoid crossing the Nyquist rate.
    pub fn generate(osc: *Oscillator, sine_table: []const f32, phase_incr: f32) f32 {
        std.debug.assert(phase_incr >= 0);
        std.debug.assert(phase_incr < 0.5);
        std.debug.assert(osc.phase >= 0 and osc.phase < 1);

        const out = sine_value(sine_table, osc.phase);

        osc.phase += phase_incr;
        while (osc.phase >= 1) osc.phase -= 1;

        return out;
    }
};

test Demodulator {
    _ = Demodulator(31, 44_100, 882);
}

test "modulate and demodulate symbols" {
    // try test_modulate_and_demodulate(31, 44_100, 882);
    // try test_modulate_and_demotulate(31, 44_100, 44.1);
    // try test_modulate_and_demotulate(31, 1_000, 1);
    // try test_modulate_and_demotulate(60, 10_000, 1);
    // try test_modulate_and_demotulate(60, 12_000, 20);
}

fn test_modulate_and_demodulate(
    comptime N: u16,
    comptime sample_rate: comptime_float,
    comptime baud: comptime_float,
) !void {
    const Mod = Modulator(N, sample_rate, baud, 64);
    const Demod = Demodulator(N, sample_rate, baud);

    var mod = Mod.init();
    var demod = Demod.init();

    const symbols = [_]Symbol{
        .waiting,            .waiting,            .waiting,            .ready,
        .{ .payload = 't' }, .{ .payload = 'a' }, .{ .payload = 'c' }, .{ .payload = 'o' },
        .{ .payload = 'h' }, .{ .payload = 'a' }, .{ .payload = 't' }, .ready,
        .ready,              .ready,              .ready,              .ready,
    };

    var signal: [Mod.symbol_len]f32 = undefined;
    var results: [symbols.len]Symbol = undefined;
    var results_index: usize = 0;

    for (symbols, 0..) |sym, i| {
        mod.modulate(sym, &signal);

        var sig: []const f32 = &signal;
        while (demod.demodulate(sig)) |res| {
            const read, const status = res;
            sig = sig[read..];
            switch (status) {
                .disconnected => {
                    std.debug.panic("could not demodulate; iter {d}, expected {}", .{ i, sym });
                },
                .symbol => |s| {
                    results[results_index] = s;
                    results_index += 1;
                },
            }
        }
    }

    // The last symbol is not expected to have been decoded because
    // the demodulator always needs some lookahead samples.
    try std.testing.expectEqual(@as(usize, 15), results_index);
    try std.testing.expectEqualSlices(Symbol, symbols[0..15], results[0..15]);
}
