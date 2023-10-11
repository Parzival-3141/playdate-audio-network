const std = @import("std");
const pdapi = @import("playdate");
const modem = @import("modem").modem2;

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
var phase: f32 = 0;
var frequency: f32 = 110;
var dc_amp: i16 = 0;

const AudioGenerator = struct {
    name: []const u8,
    func: *const fn ([*]i16, [*]i16, u32) callconv(.C) void,
};

const audio_generators = [_]AudioGenerator{
    .{ .name = "Sine", .func = generate_sine },
    .{ .name = "Sine sequence", .func = generate_sine_sequence },
    .{ .name = "Modulated data", .func = generate_modulated_data },
    .{ .name = "DC", .func = generate_dc },
};

var current_audio_generator: u4 = 0;

fn audioCallback(context: ?*anyopaque, left: [*c]i16, right: [*c]i16, len: c_int) callconv(.C) c_int {
    _ = context;

    if (sound_toggle) {
        const gen = audio_generators[current_audio_generator];
        gen.func(left, right, @intCast(len));
    } else {
        @memset(left[0..@intCast(len)], 0);
        @memset(right[0..@intCast(len)], 0);
    }

    return 1;
}

var current_amp: i16 = 0;

fn recordCallback(_: ?*anyopaque, input: [*c]i16, len: c_int) callconv(.C) c_int {
    const ulen: usize = @intCast(len);
    const samples: [*]i16 = @ptrCast(input);
    var max: i16 = 0;
    for (samples[0..ulen]) |samp| {
        var abs: i16 = samp;
        if (abs < 0) abs = -abs;
        max = @max(max, abs);
    }
    current_amp = @max(current_amp, max);
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

    if (pushed.left) current_audio_generator = @intCast((current_audio_generator -% 1) % audio_generators.len);
    if (pushed.right) current_audio_generator = @intCast((current_audio_generator +% 1) % audio_generators.len);
    if (pushed.up and sequence_periods_index < sequence_periods.len - 1) sequence_periods_index += 1;
    if (pushed.down and sequence_periods_index > 0) sequence_periods_index -= 1;

    if (pushed.b) sound_toggle = !sound_toggle;

    if (playdate.system.isCrankDocked() == 0) {
        frequency += playdate.system.getCrankChange() * 2;
        frequency = std.math.clamp(frequency, 1, 20_000);

        dc_amp = @intFromFloat((crank_angle / 180 - 1) * (std.math.maxInt(i16) - 1));
    }

    var buf = [_]u8{0} ** 128;
    var fbs = std.io.fixedBufferStream(&buf);

    fbs.writer().print(
        "accel_x: {d:.2}\naccel_y: {d:.2}\nsound_on: {}\ngenerator: {s}\n",
        .{ accel_x, accel_y, sound_toggle, audio_generators[current_audio_generator].name },
    ) catch unreachable;
    fbs.writer().print(
        "frequency: {d:.0}\nDC amp: {d}\nsequence period: {d}",
        .{ frequency, dc_amp, sequence_periods[sequence_periods_index] },
    ) catch unreachable;

    current_amp = 0;

    const len = std.mem.sliceTo(&buf, 0).len;
    _ = playdate.graphics.drawText(&buf, len, .ASCIIEncoding, 0, 0);

    // returning 1 signals to the OS to draw the frame.
    // we always want this frame drawn
    return 1;
}

fn convert_sample(x: f32) i16 {
    return @intFromFloat(x * @as(f32, @floatFromInt((std.math.maxInt(i16) - 1))));
}

fn generate_dc(left: [*]i16, right: [*]i16, count: u32) callconv(.C) void {
    for (0..count) |i| {
        left[i] = dc_amp;
        right[i] = 0;
    }
}

fn generate_sine_sample(out: *i16, incr: f32) void {
    const index_float = phase * sin_tab.len;
    const index: u32 = @intFromFloat(index_float);
    const frac = index_float - @as(f32, @floatFromInt(index));

    const y1 = sin_tab[index % sin_tab.len];
    const y2 = sin_tab[(index + 1) % sin_tab.len];
    const out_float = y1 * (1 - frac) + (y2 * frac);
    out.* = convert_sample(out_float);

    phase += incr;
    while (phase >= 1) phase -= 1;
    while (phase < 0) phase += 1;
}

fn generate_sine(left: [*]i16, right: [*]i16, count: u32) callconv(.C) void {
    for (left[0..count], right[0..count]) |*out, *dead| {
        generate_sine_sample(out, frequency / 44_100);
        dead.* = 0;
    }
}

const Modulator = modem.Modulator(31, 44_100, 882, 64);
var modulator = Modulator.init();
var modulator_sample_count: u16 = 0;
const modulator_data: []const u8 =
    \\One morning, as Gregor Samsa was waking up from anxious dreams, he
    \\discovered that in bed he had been changed into a monstrous vermin.
    \\He lay on his armour-hard back and saw, as he lifted his head up a little,
    \\his brown, arched abdomen divided up into rigid bow-like sections.
    \\From this height the blanket, just about ready to slide off completely,
    \\could hardly stay in place. His numerous legs, pitifully thin in comparis-
    \\on to the rest of his circumference, flickered helplessly before his eyes.
    \\
    \\
;
var modulator_data_index: u16 = 0;

var next_symbol_signal_buf: [Modulator.symbol_len]f32 = undefined;
var next_symbol_ready = false;
var next_symbol_index: u16 = 0;

fn generate_modulated_data(left: [*]i16, right: [*]i16, count: u32) callconv(.C) void {
    for (left[0..count], right[0..count]) |*out, *dead| {
        dead.* = 0;

        if (!next_symbol_ready) {
            var symbol: modem.Symbol = undefined;
            symbol = .{ .payload = modulator_data[modulator_data_index] };
            modulator_data_index = @intCast((modulator_data_index + 1) % modulator_data.len);

            modulator.modulate(symbol, &next_symbol_signal_buf);
            next_symbol_ready = true;
            next_symbol_index = 0;
        }

        out.* = convert_sample(next_symbol_signal_buf[next_symbol_index]);
        next_symbol_index += 1;
        if (next_symbol_index == next_symbol_signal_buf.len) {
            next_symbol_ready = false;
        }
    }
}

/// Changed via D-pad input.
var sequence_periods_index: u8 = 0;
/// How many samples a single frequency is held.
const sequence_periods = [_]u32{
    9,    18,   36,   72,    128,   256,   512,
    1024, 2756, 5512, 11025, 22050, 44100, 88200,
};

/// For each step in the sequence, multiply the base frequency by the corresponding element.
const frequency_multiplier_sequence = [_]f32{ 1, 3, 2, 4 };
/// Where we are in the sequence.
var frequency_multiplier_index: u8 = 0;
/// Total sample accumulator; when this reaches the current period len,
/// reset to zero and move to next multiplier in the sequence.
var seq_sample_count: usize = 0;

fn generate_sine_sequence(left: [*]i16, right: [*]i16, count: u32) callconv(.C) void {
    for (left[0..count], right[0..count]) |*out, *dead| {
        const mul = frequency_multiplier_sequence[frequency_multiplier_index];
        const incr = (frequency / 44_100) * mul;
        generate_sine_sample(out, incr);
        dead.* = 0;
        seq_sample_count += 1;
        if (seq_sample_count >= sequence_periods[sequence_periods_index]) {
            frequency_multiplier_index = @intCast((frequency_multiplier_index + 1) % frequency_multiplier_sequence.len);
            seq_sample_count = 0;
        }
    }
}

const sin_tab = blk: {
    var vals: [64]f32 = undefined;
    for (&vals, 0..) |*v, i| {
        v.* = @floatCast(@sin(i / @as(f64, vals.len) * std.math.pi * 2));
    }
    break :blk vals;
};
