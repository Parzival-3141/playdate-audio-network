const std = @import("std");
const dbprint = std.debug.print;
const c = @cImport({
    @cInclude("portaudio.h");
});

const demodulator = @import("modem").demodulator;

var stream: ?*c.PaStream = null;
const sample_rate: f32 = 44_100;

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
        .sampleFormat = c.paFloat32 | c.paNonInterleaved,
        .suggestedLatency = in_info.*.defaultLowInputLatency,
        .hostApiSpecificStreamInfo = null,
    };

    paAssert(c.Pa_OpenStream(
        &stream,
        &in_params,
        null,
        sample_rate,
        c.paFramesPerBufferUnspecified,
        c.paNoFlag,
        paCallback,
        null,
    ));

    paAssert(c.Pa_StartStream(stream));

    visualize_goertzel();
}

fn visualize_goertzel() void {
    while (true) {
        for (goertzel_detect_frequencies, goertzel_powers, 0..) |freq, pow, j| {
            _ = j;
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

const goertzel_N = 18;
const goertzel_analyze_period = goertzel_N;
const goertzel_detect_frequencies = blk: {
    comptime var freqs: []const f32 = &[0]f32{};
    const base_freq = sample_rate / @as(f32, goertzel_N);
    var f = base_freq;
    while (f < 20_000) : (f += base_freq) {
        freqs = freqs ++ &[1]f32{f};
    }
    break :blk freqs;
};
var sample_counter: u32 = 0;
var goertzel_buffer: [goertzel_analyze_period]f32 = undefined;
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
    _ = out_buf;
    _ = user_data;
    _ = status_flags;
    _ = time_info;

    const ins: [*]const [*]const f32 = if (in_buf) |buf| @alignCast(@ptrCast(buf)) else unreachable;

    for (0..frame_count) |i| {
        goertzel_analyze_stream(ins[0][i]);
    }

    return c.paContinue;
}

fn goertzel_analyze_stream(in_samp: f32) void {
    goertzel_buffer[sample_counter] = in_samp;
    sample_counter += 1;
    if (sample_counter == goertzel_buffer.len) {
        sample_counter = 0;

        goertzel_powers[0] = demodulator.goertzel(
            goertzel_buffer[0..],
            goertzel_detect_frequencies[0] / @as(f32, sample_rate),
            goertzel_N,
        );
        goertzel_strongest_freq_index = 0;

        for (goertzel_detect_frequencies[1..], goertzel_powers[1..], 1..) |freq, *pow, j| {
            pow.* = demodulator.goertzel(goertzel_buffer[0..], freq / @as(f32, sample_rate), goertzel_N);
            // NOTE: this can be calculated outside of audio callback,
            // but requires reading all powers atomically
            if (pow.* > goertzel_powers[j - 1]) {
                goertzel_strongest_freq_index = j;
            }
        }
    }
}

fn paAssert(code: c.PaError) void {
    if (code != c.paNoError) {
        std.debug.panic("error: portaudio: {s}", .{c.Pa_GetErrorText(code)});
    }
}
