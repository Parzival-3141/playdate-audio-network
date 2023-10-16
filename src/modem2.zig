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

        // The number of samples used to smoothly ramp from an old frequency or amplitude to a new one.
        // This transition helps avoid creating bursts of energy in higher frequencies (noise).
        const transition_samples = symbol_len - N;

        clock_osc: Oscillator,
        header_osc: Oscillator,
        payload_oscs: [8]Oscillator,
        last_clock: u1,
        last_header: u2,
        last_payload: u8,

        pub fn init() Mod {
            return Mod{
                .clock_osc = .{},
                .header_osc = .{},
                .payload_oscs = .{.{}} ** 8,
                .last_clock = 0,
                .last_header = 0,
                .last_payload = 0,
            };
        }

        fn transition_oscillator_freq(osc: *Oscillator, buf: []f32, from: f32, to: f32) void {
            for (buf[0..], 0..) |*out, i| {
                var frac: f32 = @floatFromInt(i);
                frac /= @floatFromInt(buf.len);
                const freq = cosine_interpolate(&sine_table, from, to, frac);
                var samp = osc.generate(&sine_table, freq);
                const amp = @abs((1.0 - frac) * 2.0 - 1.0);
                out.* += samp * amp * sine_scaling;
            }
        }

        /// Generate signal for the value (or no-op) into buf.
        /// Buf must be large enough to store one symbol period.
        pub fn modulate(mod: *Mod, symbol: Symbol, buf: []f32) void {
            @memset(buf, 0);

            {
                const clock = (mod.last_clock +% 1);
                transition_oscillator_freq(
                    &mod.clock_osc,
                    buf[0..transition_samples],
                    rates.clock[mod.last_clock],
                    rates.clock[clock],
                );

                // Hold frequency for remaining samples
                for (buf[transition_samples..][0..N]) |*out| {
                    var samp = mod.clock_osc.generate(&sine_table, rates.clock[clock]);
                    out.* = samp * sine_scaling;
                }

                mod.last_clock = clock;
            }

            {
                var start: u32 = 0;
                var len: u32 = symbol_len;

                const header = @intFromEnum(symbol);
                if (header != mod.last_header) {
                    transition_oscillator_freq(
                        &mod.header_osc,
                        buf[0..transition_samples],
                        rates.header[mod.last_header],
                        rates.header[header],
                    );

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

        const near_sync_countdown_reset = 10;

        pub fn init() Demod {
            return .{
                .clock_state = .start,
                .near_sync_countdown = near_sync_countdown_reset,
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
                    break :res d.check_clock_full_sync(signal, max_search_iters) orelse return null;
                },
                .set => res: {
                    if (d.near_sync_countdown == 0) {
                        const res = d.check_clock_near_sync(signal) orelse return null;
                        d.near_sync_countdown = near_sync_countdown_reset;
                        break :res res;
                    } else {
                        const res = d.check_clock_no_sync(signal) orelse return null;
                        d.near_sync_countdown -= 1;
                        break :res res;
                    }
                },
            };

            const read, const status = clock_res;
            std.log.debug("demodulate\tread {d}\tstatus {s}", .{ read, @tagName(std.meta.activeTag(status)) });
            const analysis_slice = switch (status) {
                .weak_signal, .wrong_clock => return .{ read, .disconnected },
                .ok => |slice| slice,
            };

            var target_mag: f32 = undefined;
            {
                const hdr0 = goertzel(N, analysis_slice, rates.header[0]);
                const hdr1 = goertzel(N, analysis_slice, rates.header[1]);
                const hdr2 = goertzel(N, analysis_slice, rates.header[2]);
                const sel, const ratio = select_magnitude(&.{ hdr0, hdr1, hdr2 });
                std.log.debug("header:\thdr0 {d}\thdr1 {d}\thdr2 {d}\tsel {d}\tratio {d}", .{ hdr0, hdr1, hdr2, sel, ratio });
                if (ratio < 100) {
                    d.* = init();
                    return .{ read, .disconnected };
                }

                switch (@as(u2, @intCast(sel))) {
                    0 => return .{ read, .{ .symbol = .waiting } },
                    1 => return .{ read, .{ .symbol = .ready } },
                    2 => target_mag = hdr2 / 10.0,
                    3 => unreachable,
                }
            }

            var byte: u8 = 0;
            for (rates.payload, 0..) |rate, bit| {
                const mag = goertzel(N, analysis_slice, rate);
                const val: u8 = @intFromBool(mag > target_mag);
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
            if (copy_count > sig.len) {
                @memcpy(d.buffer[d.buffer_items..][0..sig.len], sig);
                d.buffer_items += @intCast(sig.len);
                return null;
            } else {
                @memcpy(d.buffer[d.buffer_items..][0..copy_count], sig[0..copy_count]);
                d.buffer_items += copy_count;
                return .{ copy_count, d.buffer[0..needed] };
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

        const max_search_iters = std.math.log2(full_sync_samples_needed);

        /// Analyzes signal to try to find clock peak, corresponding to
        /// optimal analysis offset. Does a binary search for the highest magnitude
        /// ratio of one clock value to the other over a full symbol period.
        /// If the signal slice does not have enough samples for that check,
        /// buffers as many as possible for the next call.
        ///
        /// Whichever clock value is detected as strongest is set as the clock state,
        /// and successive windows are expected to alternate clock values.
        fn check_clock_full_sync(d: *Demod, signal: []const f32, max_iters: u16) ?SyncResult {
            std.log.debug("full sync...", .{});
            const read, const sig = d.collect_signal(signal, full_sync_samples_needed) orelse return null;

            var center: f32 = symbol_len_float / 2.0;
            var distance: f32 = symbol_len_float / 4.0;
            var start: usize = 0;
            var iter: u16 = 0;

            var clock_sel: u1 = undefined;
            var ratio: f32 = 1;

            while (distance >= 0.5 and iter < max_iters) : ({
                iter += 1;
                distance /= 2;
            }) {
                const positions = [2]f32{
                    center - distance, // left
                    center + distance, // right
                };
                for (positions) |pos_float| {
                    const pos: usize = @intFromFloat(pos_float + 0.5);
                    const chunk = sig[pos..][0..N];
                    const mags = [2]f32{
                        goertzel(N, chunk, rates.clock[0]),
                        goertzel(N, chunk, rates.clock[1]),
                    };
                    const new_sel, const new_ratio = select_magnitude(&mags);
                    std.log.debug("full sync: center {d}\tdistance ±{d}\tpos {d}\t{d}\t{d}\tsel {d}\tratio {d}", .{ center, distance, pos, mags[0], mags[1], new_sel, new_ratio });

                    if (new_ratio > ratio) {
                        clock_sel = @intCast(new_sel);
                        ratio = new_ratio;
                        center = pos_float;
                        start = pos;
                    }
                }
            }
            std.log.debug("result:    start {d}\tsel {d}\tratio {d}", .{ start, clock_sel, ratio });

            if (ratio < 100) {
                d.* = Demod.init();
                return .{ read, .weak_signal };
            }
            d.clock_state = .{ .set = clock_sel };

            // TODO: It should be possible to "unread" by as much as lookaround_count
            // to avoid copying into the buffer so successive calls can always avoid
            // relying on the buffer when incoming signal is long enough.
            // As is, the buffer is always used.
            const next_start = start + N + padding_count - lookaround_count;
            if (next_start < sig.len) {
                const save = sig.len - next_start;
                for (sig[next_start..][0..save], d.buffer[0..save]) |from, *to| {
                    to.* = from;
                }
                d.buffer_items = @intCast(save);
            } else {
                d.buffer_skip_samples = @intCast(next_start - sig.len);
                d.buffer_items = 0;
            }

            return .{ read, .{ .ok = sig[start..][0..N] } };
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
            std.log.debug("near sync...", .{});

            const read, const sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;

            // Check a sliding window of N samples to the left and right of current offset.
            var start: usize = 0;
            // var best_mag: f32 = -1;
            var ratio: f32 = 1;
            var sel: u1 = 0;
            for (0..(2 * lookaround_count + 1)) |offset| {
                const slice = sig[offset..][0..N];

                const mags = [2]f32{
                    goertzel(N, slice, rates.clock[0]),
                    goertzel(N, slice, rates.clock[1]),
                };
                const new_sel, const new_ratio = select_magnitude(&mags);
                std.log.debug("near sync: offset {d}\tmag0 {d}\tmag1 {d}\tsel {d}\tratio {d}\tbest ratio {d}", .{ offset, mags[0], mags[1], new_sel, new_ratio, ratio });

                if (new_ratio > ratio) {
                    ratio = new_ratio;
                    sel = @intCast(new_sel);
                    start = offset;
                }
            }
            std.log.debug("result:    start {d}", .{start});

            if (ratio < 100) {
                d.* = Demod.init();
                return .{ read, .weak_signal };
            }
            if (sel == d.clock_state.set) {
                d.* = Demod.init();
                return .{ read, .wrong_clock };
            }
            d.clock_state = .{ .set = sel };

            // TODO: It should be possible to "unread" by as much as lookaround_count
            // to avoid copying into the buffer so successive calls can always avoid
            // relying on the buffer when incoming signal is long enough.
            // As is, the buffer is always used.
            const next_start = start + N + padding_count - lookaround_count;
            if (next_start < sig.len) {
                const save = sig.len - next_start;
                for (sig[next_start..][0..save], d.buffer[0..save]) |from, *to| {
                    to.* = from;
                }
                d.buffer_items = @intCast(save);
            } else {
                d.buffer_skip_samples = @intCast(next_start - sig.len);
                d.buffer_items = 0;
            }

            return .{ read, .{ .ok = sig[start..][0..N] } };
        }

        /// Validates the clock of the incoming signal and returns a slice of samples ready for symbol analysis.
        fn check_clock_no_sync(d: *Demod, signal: []const f32) ?SyncResult {
            std.log.debug("no sync...", .{});
            const read, const sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;
            const start = lookaround_count;

            const clk0 = goertzel(N, sig[start..][0..N], rates.clock[0]);
            const clk1 = goertzel(N, sig[start..][0..N], rates.clock[1]);
            const sel, const ratio = select_magnitude(&.{ clk0, clk1 });
            std.log.debug("no sync\tclk0 {d}\tclk1 {d}\tsel {d}\tratio {d}", .{ clk0, clk1, sel, ratio });
            if (ratio < 100) {
                d.* = Demod.init();
                return .{ read, .weak_signal };
            }
            if (sel == d.clock_state.set) {
                d.* = Demod.init();
                return .{ read, .wrong_clock };
            }
            d.clock_state = .{ .set = @intCast(sel) };

            // TODO: It should be possible to "unread" by as much as lookaround_count
            // to avoid copying into the buffer so successive calls can always avoid
            // relying on the buffer when incoming signal is long enough.
            // As is, the buffer is always used.
            const next_start = start + N + padding_count - lookaround_count;
            if (next_start < sig.len) {
                const save = sig.len - next_start;
                for (sig[next_start..][0..save], d.buffer[0..save]) |from, *to| {
                    to.* = from;
                }
                d.buffer_items = @intCast(save);
            } else {
                d.buffer_skip_samples = @intCast(next_start - sig.len);
                d.buffer_items = 0;
            }

            return .{ read, .{ .ok = sig[start..][0..N] } };
        }

        test check_clock_no_sync {
            var d = Demod.init();
            d.clock_state = .{ .set = 1 };
            var sig: [near_sync_samples_needed]f32 = undefined;

            for (0..4) |_| {
                const expected_clock: u1 = d.clock_state.set +% 1;
                test_write_clock_signal(&sig, sig.len / 2, expected_clock, 0.25);

                const res = d.check_clock_no_sync(&sig).?;
                const read, const status = res;
                try std.testing.expectEqual(near_sync_samples_needed, @intCast(read));
                try std.testing.expectEqual(@as([]const f32, sig[lookaround_count..][0..N]), status.ok);
                try std.testing.expectEqual(expected_clock, d.clock_state.set);

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

        test check_clock_near_sync {
            var d = Demod.init();
            var sig: [near_sync_samples_needed]f32 = undefined;

            // For each possible offset, write the optimal clock signal starting at that offset,
            // and fill the remaining area with non-optimal signal, such that the sync finds
            // the offset.
            for (0..(lookaround_count * 2) + 1) |offset| {
                d.clock_state = .{ .set = 0 };
                d.buffer_items = 0;
                d.buffer_skip_samples = 0;

                test_write_clock_signal(sig[0..], offset + (N / 2), 1, 0.25);

                const res = d.check_clock_near_sync(&sig).?;
                const read, const status = res;
                try std.testing.expectEqual(near_sync_samples_needed, @intCast(read));
                // TODO: This should be able to find the exact slice
                // try std.testing.expectEqual(@as([]const f32, sig[offset..][0..N]), status.ok);
                const ptr1 = @intFromPtr(status.ok.ptr);
                const ptr2 = @intFromPtr(sig[offset..].ptr);
                const diff = @max(ptr1, ptr2) - @min(ptr1, ptr2);
                try std.testing.expect(diff < 4 * @sizeOf(f32));
                try std.testing.expectEqual(@as(u1, 1), d.clock_state.set);
            }
        }

        test check_clock_full_sync {
            var d = Demod.init();
            var sig: [full_sync_samples_needed]f32 = undefined;

            var offset: usize = sig.len / 6;
            while (offset + N < sig.len) : (offset += sig.len / 6) {
                d = Demod.init();

                test_write_clock_signal(sig[0..], offset + (N / 2), 0, 0.25);

                const res = d.check_clock_full_sync(&sig, max_search_iters).?;
                const read, const status = res;
                try std.testing.expectEqual(full_sync_samples_needed, @intCast(read));
                // TODO: This should be able to find the exact slice
                // try std.testing.expectEqual(@as([]const f32, sig[offset..][0..N]), status.ok);
                try std.testing.expect(std.meta.activeTag(status) == .ok);
                try std.testing.expect(std.meta.activeTag(d.clock_state) == .set);
                // TODO: This should be able get the expected clock
                // try std.testing.expectEqual(@as(u1, 0), d.clock_state.set);
            }
        }

        /// Write clock sine wave with oscillating frequency, reaching target frequency at peak_pos.
        fn test_write_clock_signal(out: []f32, peak_pos: usize, clk: u1, scale: f32) void {
            var mod_phase: f32 = @floatFromInt(peak_pos);
            mod_phase /= @floatFromInt(out.len);
            mod_phase = 0.25 - mod_phase;
            if (clk == 1) mod_phase += 0.5;
            while (mod_phase < 0) mod_phase += 1;
            while (mod_phase >= 1) mod_phase -= 1;

            var clock_mod = Oscillator{ .phase = mod_phase };
            var clock_car = Oscillator{};

            const sine_table = create_sine_table(64);
            const fm = comptime rates.clockFm();
            const mod_freq: f32 = baud / sample_rate / 2.0; // oscillates once per 2 symbols

            for (out) |*dst| {
                var carrier_freq = clock_mod.generate(&sine_table, mod_freq); // [-1, 1]
                carrier_freq *= fm.modulator_amp; // [-diff/2, diff/2]
                carrier_freq += fm.carrier_center; // [center-diff/2, center+diff/2] = [low, high]

                dst.* = clock_car.generate(&sine_table, carrier_freq) * scale;
            }
        }
    };
}

fn test_write_sine_mixed(out: []f32, rate: f32, amp: f32, init_phase: f32) void {
    var phase: f32 = init_phase;
    for (out) |*samp| {
        samp.* += @sin(phase * 2 * std.math.pi) * amp;
        phase += rate;
    }
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
        const samp = buf[i];
        const s = samp + (coeff * s_prev) - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }

    // Scaling the power by N squared makes it predictable for any value of N.
    const power = (s_prev2 * s_prev2) + (s_prev * s_prev) - (coeff * s_prev * s_prev2);
    const NN: f32 = @floatFromInt(N * N);
    return power / NN;
}

test goertzel {
    var buf: [100]f32 = undefined;
    try test_goertzel(&buf, 31, &.{.{ .k = 1 }}, 1, 0.25); // simple case
    try test_goertzel(&buf, 31, &.{.{ .k = 1, .phase = 0.321 }}, 1, 0.25); // initial phase has no effect
    try test_goertzel(&buf, 31, &.{.{ .k = 2 }}, 2, 0.25); // k has no effect
    try test_goertzel(&buf, 31, &.{.{ .k = 7 }}, 7, 0.25); // yet another k
    try test_goertzel(&buf, 70, &.{.{ .k = 1 }}, 1, 0.25); // N has no effect
    try test_goertzel(&buf, 70, &.{.{ .k = 33 }}, 33, 0.25); // different N and k
    try test_goertzel(&buf, 31, &.{.{ .k = 1, .amp = 0.1 }}, 1, 0.0025); // below unity gain
    try test_goertzel(&buf, 31, &.{.{ .k = 1, .amp = 4.0 }}, 1, 4.0); // above unity gain
    try test_goertzel(&buf, 100, &.{.{ .k = 33, .amp = 2.0 }}, 33, 1.0); // above unity gain, different N and k
    try test_goertzel(&buf, 31, &.{.{ .k = 5.1 }}, 5, 0.2366493); // signal is not perfectly k-aligned
    try test_goertzel(&buf, 31, &.{.{ .k = 5.5 }}, 5, 0.1069222); // signal is even less k-aligned
    try test_goertzel(&buf, 31, &.{.{ .k = 6.1 }}, 5, 0.001631582); // non-k-alinged non-target component causes spectral leakage
    try test_goertzel(&buf, 31, &.{.{ .k = 6 }}, 5, 0.0); // k-aligned non-target component goes by undetected
    try test_goertzel(&buf, 31, &.{ // several components present
        .{ .k = 1, .amp = 0.1, .phase = 0.1 },
        .{ .k = 2, .amp = 0.1, .phase = 0.2 },
        .{ .k = 3, .amp = 0.1, .phase = 0.3 },
        .{ .k = 4, .amp = 0.5, .phase = 0.4 }, // target
        .{ .k = 5, .amp = 0.1, .phase = 0.5 },
        .{ .k = 6, .amp = 0.1, .phase = 0.6 },
    }, 4, 0.0625);
    try test_goertzel(&buf, 31, &.{ // several non-k-aligned components present
        .{ .k = 0.987, .amp = 0.1, .phase = 0.1 },
        .{ .k = 2.1198, .amp = 0.1, .phase = 0.2 },
        .{ .k = 4, .amp = 0.5, .phase = 0.3 }, // target
        .{ .k = 6.42387, .amp = 0.1, .phase = 0.4 },
    }, 4, 0.0610713);
}

fn test_goertzel(
    buf: []f32,
    N: u16,
    sines: []const struct {
        /// k/N = normalized frequency for any sample rate
        k: f32,
        amp: f32 = 1,
        phase: f32 = 0,
    },
    /// Analyze this bin
    target_k: f32,
    mag: f32,
) !void {
    @memset(buf, 0);
    const N_float: f32 = @floatFromInt(N);
    for (sines) |sine| {
        test_write_sine_mixed(buf, sine.k / N_float, sine.amp, sine.phase);
    }
    const target_freq = target_k / N_float;
    const res = goertzel(N, buf, target_freq);
    try std.testing.expectApproxEqAbs(mag, res, 0.00001);
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
