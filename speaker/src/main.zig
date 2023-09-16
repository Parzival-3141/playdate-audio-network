const std = @import("std");
const dbprint = std.debug.print;
const c = @cImport({
    @cInclude("portaudio.h");
    // @cInclude("olive.c");
    @cInclude("vc.c");
});

var stream: ?*c.PaStream = null;

const WIDTH = 800;
const HEIGHT = 600;
var pixels: [WIDTH * HEIGHT]u32 = undefined;

const Color = packed struct(u32) {
    red: u8 = 0,
    green: u8 = 0,
    blue: u8 = 0,
    alpha: u8 = 0xFF,

    pub const white: Color = @bitCast(@as(u32, 0xFFFFFFFF));
};

export fn vc_render() c.Olivec_Canvas {
    var canvas = c.olivec_canvas(&pixels, WIDTH, HEIGHT, WIDTH);
    c.olivec_fill(canvas, 0);
    c.olivec_circle(
        canvas,
        WIDTH / 2,
        HEIGHT / 2,
        100,
        @bitCast(Color.white),
    );
    return canvas;
}

pub extern fn main() callconv(.C) c_int;

pub fn main2() !void {
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
        .channelCount = @intCast(1),
        .sampleFormat = c.paInt16 | c.paNonInterleaved,
        .suggestedLatency = out_info.*.defaultLowOutputLatency,
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

    std.log.info("Relaying audio from {s} to {s}...", .{ in_info.*.name, out_info.*.name });
    while (true) {}
}

var high: i16 = 0;
var low: i16 = 0;

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

    const ins: [*]const [*]const i16 = @alignCast(@ptrCast(in_buf.?));
    const outs: [*]const [*]i16 = @alignCast(@ptrCast(out_buf.?));

    for (ins[0][0..@intCast(frame_count)], outs[0][0..@intCast(frame_count)]) |in, *out| {
        if (in > high) {
            high = in;
            dbprint("high: {d} ({d:.3})\n", .{ high, @as(f32, @floatFromInt(high)) / std.math.maxInt(i16) });
        }
        if (in < low) {
            low = in;
            dbprint("low: {d} (-{d:.3})\n", .{ low, @as(f32, @floatFromInt(low)) / std.math.minInt(i16) });
        }
        out.* = in;
    }

    return c.paContinue;
}

fn paAssert(code: c.PaError) void {
    if (code != c.paNoError) {
        std.debug.panic("error: portaudio: {s}", .{c.Pa_GetErrorText(code)});
    }
}
