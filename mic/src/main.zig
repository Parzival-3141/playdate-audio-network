const std = @import("std");
const pdapi = @import("playdate_api_definitions.zig");

var g_playdate_image: *pdapi.LCDBitmap = undefined;
var playdate: *pdapi.PlaydateAPI = undefined;

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

fn audioCallback(context: ?*anyopaque, left: [*c]i16, right: [*c]i16, len: c_int) callconv(.C) c_int {
    _ = context;

    var btns: pdapi.PDButtons = undefined;
    playdate.system.getButtonState(&btns, null, null);

    const ulen: usize = @intCast(len);
    if (btns.b) {
        for (0..ulen) |i| {
            const sign = 2 * @as(i16, @intCast((@divFloor(i, 100)) % 2)) - 1;
            left[i] = sign * (std.math.maxInt(i16) - 1);
            right[i] = sign * (std.math.maxInt(i16) - 1);
        }
    } else {
        @memset(left[0..@intCast(len)], 0);
        @memset(right[0..@intCast(len)], 0);
    }

    return 1;
}

fn update_and_render(_: ?*anyopaque) callconv(.C) c_int {
    //TODO: replace with your own code!

    const to_draw = "Hello from Zig!";

    playdate.graphics.clear(@intFromEnum(pdapi.LCDSolidColor.ColorWhite));
    const pixel_width = playdate.graphics.drawText(to_draw, to_draw.len, .UTF8Encoding, pdapi.LCD_COLUMNS / 2 - 48, pdapi.LCD_ROWS / 2 + 48);
    _ = pixel_width;

    playdate.graphics.drawRotatedBitmap(
        g_playdate_image,
        pdapi.LCD_COLUMNS / 2,
        pdapi.LCD_ROWS / 2,
        playdate.system.getCrankAngle(),
        0.5,
        0.5,
        2,
        2,
    );

    var accel_x: f32 = undefined;
    var accel_y: f32 = undefined;
    playdate.system.getAccelerometer(&accel_x, &accel_y, null);

    var buf = [_]u8{0} ** 64;
    var fbs = std.io.fixedBufferStream(&buf);
    fbs.writer().print("accel_x: {d:.2}\naccel_y: {d:.2}", .{ accel_x, accel_y }) catch unreachable;

    const len = std.mem.len(@as([*:0]u8, @ptrCast(&buf)));

    _ = playdate.graphics.drawText(&buf, len, .ASCIIEncoding, 0, 0);

    // returning 1 signals to the OS to draw the frame.
    // we always want this frame drawn
    return 1;
}
