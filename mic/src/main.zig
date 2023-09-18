const std = @import("std");
const pdapi = @import("playdate_api_definitions.zig");

var g_playdate_image: *pdapi.LCDBitmap = undefined;
var playdate: *pdapi.PlaydateAPI = undefined;

pub fn panic(msg: []const u8, stack_trace: ?*std.builtin.StackTrace, _: ?usize) noreturn {
    _ = stack_trace;
    @setCold(true);

    var buf: [100]u8 = undefined;
    const len = @min(buf.len - 1, msg.len);
    @memcpy(buf[0..len], msg[0..len]);
    buf[len] = 0;
    playdate.system.@"error"(@ptrCast(&buf));
}

pub export fn eventHandler(_playdate: *pdapi.PlaydateAPI, event: pdapi.PDSystemEvent, arg: u32) callconv(.C) c_int {
    //TODO: replace with your own code!
    _ = arg;
    playdate = _playdate;

    playdate.system.logToConsole(@tagName(event));
    switch (event) {
        .EventInit => {
            g_playdate_image = playdate.graphics.loadBitmap("playdate_image", null).?;
            const font = playdate.graphics.loadFont("/System/Fonts/Asheville-Sans-14-Bold.pft", null).?;
            playdate.graphics.setFont(font);
            playdate.system.setPeripheralsEnabled(.{ .accelerometer = true });

            playdate.system.setUpdateCallback(update_and_render, playdate);
            _ = playdate.system.addMenuItem("test", menu_item, playdate);

            _ = playdate.sound.addSource(audioCallback, null, 0);
        },
        else => {},
    }
    return 0;
}

fn menu_item(_: ?*anyopaque) callconv(.C) void {
    playdate.system.logToConsole("hello menu");
}

var sound_toggle = false;

const sample_rate: f32 = 44_100;
const block_size = 9;
const base_signal_frequency = sample_rate / block_size;

var oscillators = [2]Oscillator{
    Oscillator.init(base_signal_frequency * 1, base_signal_frequency * 2),
    Oscillator.init(base_signal_frequency * 3, base_signal_frequency * 4),
};

const test_bit_pattern = [8][2]u1{
    .{ 0, 0 },
    .{ 1, 0 },
    .{ 1, 1 },
    .{ 0, 1 },
    .{ 1, 0 },
    .{ 1, 0 },
    .{ 1, 1 },
    .{ 1, 1 },
};
var bit_pattern_index: u3 = 0;
var sample_counter: u16 = 0;

fn audioCallback(context: ?*anyopaque, left: [*c]i16, right: [*c]i16, len: c_int) callconv(.C) c_int {
    _ = context;

    if (sound_toggle) {
        const ulen: u32 = @intCast(len);
        for (left, 0..ulen) |*buf, _| {
            const bits = test_bit_pattern[bit_pattern_index];

            // Mix output of each oscillator.
            var amp: f32 = 0;
            for (&oscillators, bits) |*osc, b| {
                amp += osc.generate_sample(b == 1);
            }

            const out: i16 = @intFromFloat(amp * @as(f32, @floatFromInt((std.math.maxInt(i16) - 1))));
            buf.* = out;

            // Every block_size samples, move to next bit.
            sample_counter += 1;
            if (sample_counter > block_size) {
                sample_counter = 0;
                bit_pattern_index +%= 1;
            }
        }
    } else {
        @memset(left[0..@intCast(len)], 0);
        @memset(right[0..@intCast(len)], 0);
    }

    return 1;
}

fn update_and_render(_: ?*anyopaque) callconv(.C) c_int {
    const to_draw = "Hello from Zig!";
    const crank_angle = playdate.system.getCrankAngle();

    playdate.graphics.clear(@intFromEnum(pdapi.LCDSolidColor.ColorWhite));
    const pixel_width = playdate.graphics.drawText(to_draw, to_draw.len, .UTF8Encoding, pdapi.LCD_COLUMNS / 2 - 48, pdapi.LCD_ROWS / 2 + 48);
    _ = pixel_width;

    playdate.graphics.drawRotatedBitmap(
        g_playdate_image,
        pdapi.LCD_COLUMNS / 2,
        pdapi.LCD_ROWS / 2,
        crank_angle,
        0.5,
        0.5,
        2,
        2,
    );

    var accel_x: f32 = undefined;
    var accel_y: f32 = undefined;
    playdate.system.getAccelerometer(&accel_x, &accel_y, null);

    var pushed: pdapi.PDButtons = undefined;
    playdate.system.getButtonState(null, &pushed, null);

    if (pushed.b) sound_toggle = !sound_toggle;

    var buf = [_]u8{0} ** 128;
    var fbs = std.io.fixedBufferStream(&buf);

    fbs.writer().print(
        "accel_x: {d:.2}\naccel_y: {d:.2}\nsound_on: {}\n",
        .{ accel_x, accel_y, sound_toggle },
    ) catch unreachable;

    const len = std.mem.sliceTo(&buf, 0).len;
    _ = playdate.graphics.drawText(&buf, len, .ASCIIEncoding, 0, 0);

    // returning 1 signals to the OS to draw the frame.
    // we always want this frame drawn
    return 1;
}

/// Stateful sine wave oscillator.
const Oscillator = struct {
    /// The amount to increase the phase on each sample.
    /// Equals frequency / sample_rate.
    /// Index 0 is increment for low frequency, index 1 for high.
    phase_incr: [2]f32,
    /// Current phase; changes every sample.
    phase: f32,

    fn init(low_freq: f32, high_freq: f32) Oscillator {
        if (low_freq <= 0 or high_freq <= 0) {
            @panic("frequencies must be > 0");
        }
        return .{
            .phase_incr = .{ low_freq / sample_rate, high_freq / sample_rate },
            .phase = 0,
        };
    }

    /// Calculates the next amplitude value and updates phase.
    /// If high is true, generates at high frequency, otherwise at low frequency.
    fn generate_sample(osc: *Oscillator, high: bool) f32 {
        const index_float = osc.phase * sin_tab.len;
        const index: u8 = @intFromFloat(index_float);
        const frac = index_float - @as(f32, @floatFromInt(index));

        const y1 = sin_tab[index % sin_tab.len];
        const y2 = sin_tab[(index + 1) % sin_tab.len];
        const out = y1 * (1 - frac) + (y2 * frac);

        osc.phase += osc.phase_incr[@intFromBool(high)];
        while (osc.phase >= 1) osc.phase -= 1;

        return out;
    }
};

const sin_tab = blk: {
    var vals: [64]f32 = undefined;
    for (&vals, 0..) |*v, i| {
        v.* = @floatCast(@sin(i / @as(f64, vals.len) * std.math.pi * 2));
        // Attenuate source lookup table so multiple oscillators
        // won't exceed 1.0 when mixed.
        v.* /= oscillators.len;
    }
    break :blk vals;
};
