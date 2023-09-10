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

    const in_params = c.PaStreamParameters{
        .device = in_device,
        .channelCount = @intCast(1),
        .sampleFormat = c.paFloat32 | c.paNonInterleaved,
        .suggestedLatency = c.Pa_GetDeviceInfo(in_device).*.defaultLowInputLatency,
        .hostApiSpecificStreamInfo = null,
    };

    const out_device = c.Pa_GetDefaultOutputDevice();
    if (out_device == c.paNoDevice) {
        @panic("No default output device!\n");
    }

    const out_params = c.PaStreamParameters{
        .device = out_device,
        .channelCount = @intCast(1),
        .sampleFormat = c.paFloat32 | c.paNonInterleaved,
        .suggestedLatency = c.Pa_GetDeviceInfo(out_device).*.defaultLowOutputLatency,
        .hostApiSpecificStreamInfo = null,
    };

    dbprint("are in/out the same: {any}\n", .{in_device == out_device});

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

    std.log.info("listening for audio..", .{});
    while (true) {}
}

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

    const ins: [*]const [*]const f32 = if (in_buf) |buf| @alignCast(@ptrCast(buf)) else unreachable;
    const outs: [*]const [*]f32 = if (out_buf) |buf| @alignCast(@ptrCast(buf)) else unreachable;

    @memcpy(outs[0][0..@intCast(frame_count)], ins[0][0..@intCast(frame_count)]);

    // var avg_sample_value: f32 = 0;
    // for (ins[0][0..@intCast(frame_count)]) |sample| {
    //     avg_sample_value += sample;
    // }

    // avg_sample_value /= @floatFromInt(frame_count);

    // dbprint("avg. sample value: {d:.3}\n", .{avg_sample_value});

    return c.paContinue;
}

fn paAssert(code: c.PaError) void {
    if (code != c.paNoError) {
        std.debug.panic("error: portaudio: {s}", .{c.Pa_GetErrorText(code)});
    }
}
