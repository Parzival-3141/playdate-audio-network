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
        .sampleFormat = c.paFloat32 | c.paNonInterleaved,
        .suggestedLatency = in_info.*.defaultHighInputLatency,
        .hostApiSpecificStreamInfo = null,
    };

    const out_device = c.Pa_GetDefaultOutputDevice();
    if (out_device == c.paNoDevice) {
        @panic("No default output device!\n");
    }
    // const out_info = c.Pa_GetDeviceInfo(out_device);
    // const out_params = c.PaStreamParameters{
    //     .device = out_device,
    //     .channelCount = @intCast(1),
    //     .sampleFormat = c.paFloat32 | c.paNonInterleaved,
    //     .suggestedLatency = out_info.*.defaultLowOutputLatency,
    //     .hostApiSpecificStreamInfo = null,
    // };

    paAssert(c.Pa_OpenStream(
        &stream,
        &in_params,
        // &out_params,
        null,
        44_100,
        c.paFramesPerBufferUnspecified,
        c.paNoFlag,
        paCallback,
        null,
    ));

    paAssert(c.Pa_StartStream(stream));

    // std.log.info("Relaying audio from {s} to {s}...", .{ in_info.*.name, out_info.*.name });
    std.log.info("Reading audio from {s}...", .{in_info.*.name});
    while (true) {}
}

// https://en.wikipedia.org/wiki/Goertzel_algorithm
// https://www.embedded.com/the-goertzel-algorithm/

/// Number of samples per block
const block_size = 9;
/// The number of frequencies the algo will operate on
const freq_components = 4;

fn goertzel() [freq_components]f32 {
    var results = [_]f32{0} ** freq_components;

    for (1..freq_components + 1) |k| {
        const w = (2 * std.math.pi / @as(f32, freq_components)) * @as(f32, @floatFromInt(k));
        const cos = @cos(w);
        const sin = @sin(w);
        const coeff = 2 * cos;

        var current: f32 = 0;
        var previous: f32 = 0;
        var prev_previous: f32 = 0;
        for (sample_bin) |sample| {
            current = coeff * previous - prev_previous + sample;
            prev_previous = previous;
            previous = current;
        }

        const real = previous - prev_previous * cos;
        const imaginary = prev_previous * sin;
        const sqr_magnitude = real * real + imaginary * imaginary;

        results[k - 1] = sqr_magnitude;
    }

    return results;
}

var sample_bin: [block_size]f32 = undefined;
var bin_index: u32 = 0;

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

    const ins: [*]const [*]const f32 = @alignCast(@ptrCast(in_buf.?));
    // const outs: [*]const [*]f32 = @alignCast(@ptrCast(out_buf.?));

    for (ins[0][0..frame_count]) |sample| {
        sample_bin[bin_index] = sample;
        bin_index += 1;
        if (bin_index == sample_bin.len) {
            const sqr_magnitudes = goertzel();

            dbprint(
                "{d:.4} | {d:.4} | {d:.4} | {d:.4}\n",
                .{ @round(sqr_magnitudes[0]), @round(sqr_magnitudes[1]), @round(sqr_magnitudes[2]), @round(sqr_magnitudes[3]) },
            );

            bin_index = 0;
        }
    }

    return c.paContinue;
}

fn paAssert(code: c.PaError) void {
    if (code != c.paNoError) {
        std.debug.panic("error: portaudio: {s}", .{c.Pa_GetErrorText(code)});
    }
}
