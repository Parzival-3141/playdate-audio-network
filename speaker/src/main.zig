const std = @import("std");
const dbprint = std.debug.print;
const c = @cImport({
    @cInclude("portaudio.h");
});

var stream: ?*c.PaStream = null;

pub fn main() !void {
    paAssert(c.Pa_Initialize());
    defer paAssert(c.Pa_Terminate());

    const in_device = c.Pa_GetDefaultInputDevice();
    if (in_device == c.paNoDevice) {
        @panic("No default input device!\n");
    }
    const in_info = c.Pa_GetDeviceInfo(in_device);
    const in_params = c.PaStreamParameters{
        .device = in_device,
        .channelCount = @intCast(1),
        .sampleFormat = c.paInt16 | c.paNonInterleaved,
        .suggestedLatency = in_info.*.defaultLowInputLatency,
        .hostApiSpecificStreamInfo = null,
    };

    const out_device = c.Pa_GetDefaultOutputDevice();
    if (out_device == c.paNoDevice) {
        @panic("No default output device!\n");
    }
    const out_info = c.Pa_GetDeviceInfo(out_device);
    const out_params = c.PaStreamParameters{
        .device = out_device,
        .channelCount = @intCast(2),
        .sampleFormat = c.paInt16 | c.paNonInterleaved,
        .suggestedLatency = out_info.*.defaultLowOutputLatency,
        .hostApiSpecificStreamInfo = null,
    };

    paAssert(c.Pa_OpenStream(
        &stream,
        &in_params,
        &out_params,
        44_100,
        c.paFramesPerBufferUnspecified,
        c.paNoFlag,
        paCallback,
        null,
    ));

    paAssert(c.Pa_StartStream(stream));

    while (true) {
        for (goertzel_detect_frequencies, goertzel_powers) |freq, pow| {
            const max_stars = 40;
            var star_buf = [_]u8{' '} ** max_stars;
            const stars = @min(@as(u32, @intFromFloat(pow * @as(f32, max_stars))), max_stars);
            for (0..stars) |i| {
                star_buf[i] = '*';
            }
            std.debug.print("{d: <8.1} Hz: {d: >10.5} {s}\n", .{ freq, pow, star_buf[0..] });
        }
        std.debug.print("strongest freq: {d:.1}\n", .{goertzel_detect_frequencies[goertzel_strongest_freq_index]});

        std.time.sleep(100 * std.time.ns_per_ms);

        for (0..goertzel_detect_frequencies.len + 1) |_| {
            // Clear line, cursor up one line
            std.debug.print("\x1b[2K\x1b[1A", .{});
        }
        // Final carriage return to return to first column
        std.debug.print("\r", .{});
    }
}

var high: i16 = 0;
var low: i16 = 0;

const sample_rate = 44_100;
const goertzel_N = 9;
const goertzel_detect_frequencies = [_]f32{
    sample_rate / @as(f32, goertzel_N) * 1,
    sample_rate / @as(f32, goertzel_N) * 2,
    sample_rate / @as(f32, goertzel_N) * 3,
    sample_rate / @as(f32, goertzel_N) * 4,
};
var sample_counter: u32 = 0;
var goertzel_buffer: [goertzel_N]f32 = undefined;
var goertzel_powers = [1]f32{0} ** goertzel_detect_frequencies.len;
var goertzel_strongest_freq_index: usize = 0;

fn paCallback(
    in_buf: ?*const anyopaque,
    out_buf: ?*anyopaque,
    frame_count: c_ulong,
    time_info: [*c]const c.PaStreamCallbackTimeInfo,
    status_flags: c.PaStreamCallbackFlags,
    user_data: ?*anyopaque,
) callconv(.C) c_int {
    _ = user_data;
    _ = status_flags;
    _ = time_info;

    const ins: [*]const [*]const i16 = if (in_buf) |buf| @alignCast(@ptrCast(buf)) else unreachable;
    const outs: [*]const [*]i16 = if (out_buf) |buf| @alignCast(@ptrCast(buf)) else unreachable;

    // Silence output
    @memset(outs[0][0..frame_count], 0);
    @memset(outs[1][0..frame_count], 0);

    for (0..frame_count) |i| {
        goertzel_buffer[sample_counter] = convertSample(ins[0][i]);
        sample_counter += 1;
        if (sample_counter == goertzel_N) {
            sample_counter = 0;
            for (goertzel_detect_frequencies, &goertzel_powers, 0..) |freq, *pow, j| {
                pow.* = goertzel(goertzel_buffer[0..], freq, goertzel_N);
                if (j == 0) {
                    goertzel_strongest_freq_index = 0;
                } else {
                    // NOTE: this can be calculated outside of audio callback,
                    // but requires reading all powers atomically
                    if (pow.* > goertzel_powers[j - 1]) {
                        goertzel_strongest_freq_index = j;
                    }
                }
            }
        }
    }

    return c.paContinue;
}

fn convertSample(val: i16) f32 {
    return @as(f32, @floatFromInt(val)) / @as(f32, @floatFromInt(std.math.maxInt(i16) - 1));
}

/// Calculates the power of the given frequency within the input sample slice.
/// N is the number of DFT terms; frequency detection is strongest when freq is
/// an integer multiple of (1/N * sample_rate).
///
/// Reference: https://en.wikipedia.org/wiki/Goertzel_algorithm#Power-spectrum_terms
fn goertzel(in: []const f32, freq: f32, N: u16) f32 {
    const freq_normalized = freq / @as(f32, sample_rate);
    const coeff = 2 * @cos(2.0 * std.math.pi * freq_normalized);

    var s_prev: f32 = 0;
    var s_prev2: f32 = 0;

    for (0..N) |i| {
        const s = in[i] + (coeff * s_prev) - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }

    const power = (s_prev2 * s_prev2) + (s_prev * s_prev) - (coeff * s_prev * s_prev2);
    // TODO: scale down based on size of N
    return power;
}

const sin_tab = blk: {
    var vals: [256]f32 = undefined;
    for (&vals, 0..) |*v, i| {
        v.* = @floatCast(@sin(i / @as(f64, vals.len) * std.math.pi * 2));
    }
    break :blk vals;
};

var phase: f32 = 0;
const frequency: f32 = 19_294;
var phase_incr: f32 = frequency / 44_100.0;

fn paAssert(code: c.PaError) void {
    if (code != c.paNoError) {
        std.debug.panic("error: portaudio: {s}", .{c.Pa_GetErrorText(code)});
    }
}
