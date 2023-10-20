const std = @import("std");
const assert = std.debug.assert;

const modem = @import("modem.zig");
const OscillatorRates = modem.OscillatorRates;
const Symbol = modem.Symbol;
const goertzel = modem.goertzel;

pub fn Demodulator(
    comptime N: u16,
    comptime sample_rate: comptime_float,
    comptime baud: comptime_float,
) type {
    comptime {
        assert(N > 0);
        assert(sample_rate > 0);
        assert(baud > 0);
    }

    const rates = OscillatorRates.init(N, sample_rate);

    return struct {
        const Demod = @This();

        const symbol_len_float = @floor(sample_rate / baud);
        pub const symbol_len: u16 = @intFromFloat(symbol_len_float);
        comptime {
            assert(symbol_len == symbol_len_float);
            if (N > symbol_len) {
                @compileError("symbol period must be greater than N samples; choose a lower baud");
            }
        }

        const lookaround_count = 2;
        const near_sync_samples_needed = N + (2 * lookaround_count);
        const full_sync_samples_needed = symbol_len + N;
        comptime {
            assert(lookaround_count < N);
            assert(near_sync_samples_needed < symbol_len);
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

                if (ratio < 50) {
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
            }

            return .{ read, .{ .symbol = .{ .payload = byte } } };
        }

        const SyncResult = struct {
            usize,
            []const f32,
        };

        fn get_slice_full_sync(d: *Demod, signal: []const f32) ?SyncResult {
            const read, const sig = d.collect_signal(signal, full_sync_samples_needed) orelse return null;
            const start = get_best_power(sig, N, N / 3); // TODO: consider using binary search
            d.buffer_finish(sig, start + symbol_len - lookaround_count);

            return .{ read, sig[start..][0..N] };
        }

        fn get_slice_near_sync(d: *Demod, signal: []const f32) ?SyncResult {
            const read, const sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;
            const start = get_best_power(sig, N, 1);
            d.buffer_finish(sig, start + symbol_len - lookaround_count);

            return .{ read, sig[start..][0..N] };
        }

        fn get_slice_no_sync(d: *Demod, signal: []const f32) ?SyncResult {
            const read, const sig = d.collect_signal(signal, near_sync_samples_needed) orelse return null;
            const start = lookaround_count;
            d.buffer_finish(sig, start + symbol_len - lookaround_count);

            return .{ read, sig[start..][0..N] };
        }

        /// When successful, returns the amount read and a slice of `needed` samples.
        /// If signal is long enough and internal buffering is not already in progress,
        /// the signal slice is used directly to avoid copying.
        /// Otherwise, buffers incoming signal until `needed` samples has been collected.
        /// Until enough samples are ready, returns null, implying `read` is equal to signal.len.
        fn collect_signal(d: *Demod, signal: []const f32, needed: usize) ?struct { usize, []const f32 } {
            if (d.buffer_skip_samples >= signal.len) {
                d.buffer_skip_samples -= @intCast(signal.len);
                return null;
            }

            const skipped: usize = d.buffer_skip_samples;
            const sig = signal[d.buffer_skip_samples..];
            d.buffer_skip_samples = 0;

            // Happy path: no buffering was needed.
            if (d.buffer_items == 0 and sig.len >= needed) {
                return .{ skipped + needed, sig[0..needed] };
            }

            // Otherwise, either the signal slice is too short or we need to append to
            // an in-progress buffering.
            const remaining: u16 = @intCast(needed - d.buffer_items);
            if (remaining > sig.len) {
                @memcpy(d.buffer[d.buffer_items..][0..sig.len], sig);
                d.buffer_items += @intCast(sig.len);
                return null;
            } else {
                @memcpy(d.buffer[d.buffer_items..][0..remaining], sig[0..remaining]);
                d.buffer_items += remaining;
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

        fn buffer_finish(d: *Demod, sig: []const f32, next_start: u16) void {
            assert(d.buffer_skip_samples == 0);

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
    assert(incr > 0);

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

    return best_pos;
}

/// Returns the index of the highest value and the ratio of it to the second highest value.
/// Compared values are clamped to a reasonable range to prevent NaN and Inf.
fn select_magnitude(mags: []const f32) struct { u16, f32 } {
    assert(mags.len > 0);

    var sel: u16 = 0;
    var max_mag: f32 = 0;
    var next_mag: f32 = 0;

    for (mags, 0..) |mag, i| {
        assert(mag >= 0);
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

test Demodulator {
    _ = Demodulator(31, 44_100, 882);
}
