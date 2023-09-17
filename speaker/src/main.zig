const std = @import("std");
const dbprint = std.debug.print;
const c = @cImport({
    @cInclude("portaudio.h");
});

var stream: ?*c.PaStream = null;

pub fn main() !void {
    paAssert(c.Pa_Initialize());
    defer paAssert(c.Pa_Terminate());

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
        null,
        &out_params,
        44_100,
        c.paFramesPerBufferUnspecified,
        c.paNoFlag,
        paCallback,
        null,
    ));

    paAssert(c.Pa_StartStream(stream));

    std.log.info("Sending signal to {s}...", .{out_info.*.name});
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
    _ = in_buf;
    _ = user_data;
    _ = status_flags;
    _ = time_info;

    const outs: [*]const [*]i16 = if (out_buf) |buf| @alignCast(@ptrCast(buf)) else unreachable;

    for (0..frame_count) |i| {
        const index_float = phase * sin_tab.len;
        const index: u32 = @intFromFloat(index_float);
        const frac = index_float - @as(f32, @floatFromInt(index));

        const y1 = sin_tab[index % sin_tab.len];
        const y2 = sin_tab[(index + 1) % sin_tab.len];
        const out_float = y1 * (1 - frac) + (y2 * frac);
        const out: i16 = @intFromFloat(out_float * @as(f32, @floatFromInt((std.math.maxInt(i16) - 1))));

        outs[0][i] = out;
        outs[1][i] = out;

        phase += phase_incr;
        while (phase >= 1) phase -= 1;
        while (phase < 0) phase += 1;
    }

    return c.paContinue;
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
