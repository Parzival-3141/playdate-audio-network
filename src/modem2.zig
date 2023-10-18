const std = @import("std");

pub const OscillatorRates = extern struct {
    header: [3]f32,
    payload: [8]f32,

    const osc_count = 11;

    /// Returns a set of oscillator phase rates (inverse frequencies) that are optimal for
    /// Goertzel detection with N DFT terms at the given sample rate.
    ///
    /// One symbol requires 11 unique frequencies, in increasing order:
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

/// Values that one side (client) can send to the other (server).
/// Both sides may simultaneously act as client and server to one another.
pub const Symbol = union(enum) {
    /// Client is running but has not detected a valid signal from the server.
    /// Server should not send a payload back while the client is in this state.
    waiting,
    /// Client has detected the server but there is no data to send.
    ready,
    /// Client has detected the server and is sending a byte of data.
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
    comptime {
        std.debug.assert(N > 0);
        std.debug.assert(sample_rate > 0);
        std.debug.assert(baud > 0);
    }

    const rates = OscillatorRates.init(N, sample_rate);

    return struct {
        const Mod = @This();

        const symbol_len_float = @floor(sample_rate / baud);
        pub const symbol_len: u16 = @intFromFloat(symbol_len_float);
        comptime {
            std.debug.assert(symbol_len == symbol_len_float);
            if (N > symbol_len) {
                @compileError("symbol period must be greater than N samples");
            }
        }

        pub const sine_table = create_sine_table(sine_table_len);

        header_osc: Oscillator,
        payload_oscs: [8]Oscillator,

        pub fn init() Mod {
            return Mod{
                .header_osc = .{ .phase = 0 },
                .payload_oscs = .{
                    .{ .phase = 0.1 },
                    .{ .phase = 0.2 },
                    .{ .phase = 0.3 },
                    .{ .phase = 0.4 },
                    .{ .phase = 0.5 },
                    .{ .phase = 0.6 },
                    .{ .phase = 0.7 },
                    .{ .phase = 0.8 },
                },
            };
        }

        /// Generate signal for the symbol into buf.
        /// Buf must be large enough to store one symbol period.
        pub fn modulate(mod: *Mod, symbol: Symbol, buf: []f32) void {
            const header = @intFromEnum(symbol);
            for (buf[0..symbol_len]) |*out| {
                var samp = mod.header_osc.generate(&sine_table, rates.header[header]);
                out.* = samp;
            }

            const byte = switch (symbol) {
                .payload => |b| b,
                else => 0,
            };

            for (0..8) |bit| {
                const mask: u8 = @as(u8, 1) << @intCast(bit);
                const bit_set = byte & mask != 0;
                if (!bit_set) {
                    continue;
                }

                for (buf[0..symbol_len]) |*out| {
                    var samp = mod.payload_oscs[bit].generate(&sine_table, rates.payload[bit]);
                    out.* += samp;
                }
            }

            // Total scale factor is inverse of max number of
            // simultaneously mixed oscillators to avoid clipping
            const scale: f32 = 1.0 / OscillatorRates.osc_count;

            // Shape symbol window with a full sine cycle.
            for (buf[0..symbol_len], 0..) |*out, i| {
                var frac = @as(f32, @floatFromInt(i)) / symbol_len_float + 0.75;
                if (frac >= 1) frac -= 1;
                var amp = sine_value(&sine_table, frac); // [-1, 1]
                amp *= 0.5; // [-0.5, 0.5]
                amp += 0.5; // [0, 1]
                amp *= scale;
                out.* *= amp;
            }
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

pub fn Demodulator(
    comptime N: u16,
    comptime sample_rate: comptime_float,
    comptime baud: comptime_float,
) type {
    comptime {
        std.debug.assert(N > 0);
        std.debug.assert(sample_rate > 0);
        std.debug.assert(baud > 0);
    }

    const rates = OscillatorRates.init(N, sample_rate);

    return struct {
        const Demod = @This();

        const symbol_len_float = @floor(sample_rate / baud);
        pub const symbol_len: u16 = @intFromFloat(symbol_len_float);
        comptime {
            std.debug.assert(symbol_len == symbol_len_float);
            if (N > symbol_len) {
                @compileError("symbol period must be greater than N samples; choose a lower baud");
            }
        }

        const lookaround_count = 2;
        const near_sync_samples_needed = N + (2 * lookaround_count);
        const full_sync_samples_needed = symbol_len + N;
        comptime {
            std.debug.assert(lookaround_count < N);
            std.debug.assert(near_sync_samples_needed < symbol_len);
        }

        const near_sync_countdown_reset = 6;

        pub const SyncState = union(enum) {
            start,
            near_sync_countdown: u16,
        };

        sync_state: SyncState,
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
                .sync_state = .start,
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
                /// The incoming signal did not yield a clear signal.
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
        /// the end of `signal` was reached, i.e. `read` is implicitly equal to signal.len.
        pub fn demodulate(d: *Demod, signal: []const f32) ?Result {
            std.log.debug("demodulating symbol", .{});

            const read, const analysis_slice = switch (d.sync_state) {
                .start => res: {
                    const res = d.get_slice_full_sync(signal) orelse return null;
                    d.sync_state = .{ .near_sync_countdown = near_sync_countdown_reset };
                    break :res res;
                },
                .near_sync_countdown => |count| res: {
                    if (count == 0) {
                        const res = d.get_slice_near_sync(signal) orelse return null;
                        d.sync_state.near_sync_countdown = near_sync_countdown_reset;
                        break :res res;
                    } else {
                        const res = d.get_slice_no_sync(signal) orelse return null;
                        d.sync_state.near_sync_countdown -= 1;
                        break :res res;
                    }
                },
            };

            const target_mag = blk: {
                const hdr0 = goertzel(N, analysis_slice, rates.header[0]);
                const hdr1 = goertzel(N, analysis_slice, rates.header[1]);
                const hdr2 = goertzel(N, analysis_slice, rates.header[2]);
                const sel, const ratio = select_magnitude(&.{ hdr0, hdr1, hdr2 });
                std.log.debug("header:\thdr0 {d}\thdr1 {d}\thdr2 {d}\tsel {d}\tratio {d}", .{ hdr0, hdr1, hdr2, sel, ratio });
                if (ratio < 50) {
                    std.log.debug("{d}", .{analysis_slice});
                    d.* = init();
                    return .{ read, .disconnected };
                }

                switch (@as(u2, @intCast(sel))) {
                    0 => return .{ read, .{ .symbol = .waiting } },
                    1 => return .{ read, .{ .symbol = .ready } },
                    2 => break :blk hdr2 / 10.0,
                    3 => unreachable,
                }
            };

            var byte: u8 = 0;
            for (rates.payload, 0..) |rate, bit| {
                const mag = goertzel(N, analysis_slice, rate);
                const val: u8 = @intFromBool(mag > target_mag);
                const mask: u8 = val << @intCast(bit);
                byte |= mask;

                std.log.debug("bit {d} mag {d:0>.10} target {d:0>.10}{s}", .{ bit, mag, target_mag, if (val == 0) "" else " *" });
            }

            return .{ read, .{ .symbol = .{ .payload = byte } } };
        }

        const SyncResult = struct {
            usize,
            []const f32,
        };

        fn get_slice_full_sync(d: *Demod, signal: []const f32) ?SyncResult {
            std.log.debug("full sync...", .{});

            const read, const sig = d.collect_signal(signal, full_sync_samples_needed) orelse return null;
            const start = get_best_power(sig, N, N / 3); // TODO: consider using binary search
            std.log.debug("best position: {}", .{start});
            d.buffer_finish(sig, start, start + symbol_len - lookaround_count);

            return .{ read, sig[start..][0..N] };
        }

        fn get_slice_near_sync(d: *Demod, signal: []const f32) ?SyncResult {
            std.log.debug("near sync...", .{});

            const read, const sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;
            const start = get_best_power(sig, N, 1);
            std.log.debug("best position: {}", .{start});
            d.buffer_finish(sig, start, start + symbol_len - lookaround_count);

            return .{ read, sig[start..][0..N] };
        }

        fn get_slice_no_sync(d: *Demod, signal: []const f32) ?SyncResult {
            std.log.debug("no sync...", .{});

            const read, const sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;
            const start = lookaround_count;
            d.buffer_finish(sig, start, start + symbol_len - lookaround_count);

            return .{ read, sig[start..][0..N] };
        }

        fn write_sig_to_stdout(sig: []const f32) void {
            for (sig) |float| {
                const buf_bytes: [4]u8 = @bitCast(float);
                std.io.getStdOut().writeAll(&buf_bytes) catch unreachable;
            }
        }

        /// When successful, returns the amount read and a slice of `needed` samples.
        /// If signal is long enough and internal buffering is not already in progress,
        /// the signal slice is used directly to avoid copying.
        /// Otherwise, buffers incoming signal until `needed` samples has been collected.
        /// Until enough samples are ready, returns null, implying `read` is equal to signal.len.
        fn collect_signal(d: *Demod, signal: []const f32, needed: usize) ?struct { usize, []const f32 } {
            if (d.buffer_skip_samples >= signal.len) {
                // write_sig_to_stdout(signal); // TEMP
                // write_sig_to_stdout(&.{ -1, 0, -1, 0, -1, 0, -1, 0 }); // TEMP
                std.log.debug("collect: skipped {} (entire slice)", .{signal.len});

                d.buffer_skip_samples -= @intCast(signal.len);
                return null;
            }

            const skipped: usize = d.buffer_skip_samples;
            const sig = signal[d.buffer_skip_samples..];
            d.buffer_skip_samples = 0;
            if (skipped > 0) {
                // write_sig_to_stdout(signal[0..skipped]); // TEMP
                // write_sig_to_stdout(&.{ 1, 0, 1, 0, 1, 0, 1, 0 }); // TEMP
                std.log.debug("collect: skipped {}, signal remainder {}", .{ skipped, sig.len });
            }

            // Happy path: no buffering was needed.
            if (d.buffer_items == 0 and sig.len >= needed) {
                std.log.debug("collect: signal had all needed ({})", .{needed});
                // write_sig_to_stdout(sig[0..needed]); // TEMP
                // write_sig_to_stdout(&.{ 1, 1, 1, 1, 1, 1, 1, 1 }); // TEMP

                return .{ skipped + needed, sig[0..needed] };
            }

            // Otherwise, either the signal slice is too short or we need to append to
            // an in-progress buffering.
            const remaining: u16 = @intCast(needed - d.buffer_items);
            if (remaining > sig.len) {
                @memcpy(d.buffer[d.buffer_items..][0..sig.len], sig);
                d.buffer_items += @intCast(sig.len);

                std.log.debug("collect: still need {}, copying {} into buffer", .{ remaining, sig.len });
                // write_sig_to_stdout(signal); // TEMP
                // write_sig_to_stdout(&.{ -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75 }); // TEMP

                return null;
            } else {
                @memcpy(d.buffer[d.buffer_items..][0..remaining], sig[0..remaining]);
                d.buffer_items += remaining;

                std.log.debug("collect: copying remaining {} into buffer and returning", .{remaining});
                // write_sig_to_stdout(signal[0..remaining]); // TEMP
                // write_sig_to_stdout(&.{ 1, 0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75 }); // TEMP

                return .{ skipped + remaining, d.buffer[0..needed] };
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

        fn buffer_finish(d: *Demod, sig: []const f32, start: u16, next_start: u16) void {
            std.debug.assert(d.buffer_skip_samples == 0);
            std.log.debug("prepare buffer: sig len {} start {} next {}", .{ sig.len, start, next_start });
            defer std.log.debug("buffer has items {} skip={}", .{ d.buffer_items, d.buffer_skip_samples });

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
        }
    };
}

fn get_best_power(sig: []const f32, N: u16, incr: u16) u16 {
    std.debug.assert(incr > 0);

    var pos: u16 = 0;
    var best_pos = pos;
    var best_power: f32 = 0;
    while (pos + N < sig.len) : (pos += incr) {
        var power: f32 = 0;
        for (sig[pos..][0..N]) |samp| {
            power += samp * samp;
        }
        if (power > best_power) {
            best_pos = pos;
            best_power = power;
        }
    }
    std.log.debug("best power {d}", .{best_power});
    return best_pos;
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

    next_mag = std.math.clamp(next_mag, 0.000001, 1000.0);
    max_mag = std.math.clamp(max_mag, 0.000001, 1000.0);

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
        var phase: f32 = sine.phase;
        for (buf) |*samp| {
            samp.* += @sin(phase * 2 * std.math.pi) * sine.amp;
            phase += sine.k / N_float;
        }
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
    try test_modulate_and_demodulate(26, 1_000, 10);
    try test_modulate_and_demodulate(26, 44_100, 100);
    try test_modulate_and_demodulate(26, 44_100, 441);
    try test_modulate_and_demodulate(30, 48_000, 500);
    try test_modulate_and_demodulate(40, 10_000, 50);
    try test_modulate_and_demodulate(40, 14_000, 160);
    try test_modulate_and_demodulate(60, 48_000, 10);
    try test_modulate_and_demodulate(60, 48_000, 50);
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

        var pos: usize = 0;
        while (pos < signal.len) {
            const read, const status = demod.demodulate(signal[pos..]) orelse break;
            pos += read;
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

    try std.testing.expectEqual(symbols.len, results_index);
    try std.testing.expectEqualSlices(Symbol, &symbols, &results);
}
