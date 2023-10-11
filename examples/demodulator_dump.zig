const std = @import("std");
const dbprint = std.debug.print;
const c = @cImport({
    @cInclude("portaudio.h");
});

const options = @import("options");
pub const std_options = struct {
    pub const log_level = switch (options.log_level) {
        .err => .err,
        .warn => .warn,
        .info => .info,
        .debug => .debug,
    };
};

const modem = @import("modem").modem2;

var stream: ?*c.PaStream = null;
const sample_rate: f32 = 44_100;

pub fn main() !void {
    paAssert(c.Pa_Initialize());
    defer paAssert(c.Pa_Terminate());

    // const in_device = c.Pa_GetDefaultInputDevice();
    const in_device = try get_input_device();
    if (in_device == c.paNoDevice) {
        @panic("invalid input device!\n");
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

    dump_bitstream();
}

fn get_input_device() !c.PaDeviceIndex {
    const num_devices: usize = num: {
        const count = c.Pa_GetDeviceCount();
        paAssert(@min(count, 0)); // portaudio errors are negative values, so we just check for that
        break :num @intCast(count);
    };
    if (num_devices == 0) return error.NoAudioDevices;

    const default = d: {
        const dev = c.Pa_GetDefaultInputDevice();
        break :d if (dev != c.paNoDevice) dev else null;
    };

    const stdout = std.io.getStdOut().writer();

    try stdout.writeAll("\nChoose an input device to listen on (will use default if none is selected):\n");

    if (default) |d| {
        const info: *const c.PaDeviceInfo = c.Pa_GetDeviceInfo(d);
        try stdout.print("  default: {s}\n", .{info.name});
    }

    var i: usize = 1;
    for (0..num_devices) |d| {
        const info: *const c.PaDeviceInfo = c.Pa_GetDeviceInfo(@intCast(d));
        try stdout.print("  {d}: {s}\n", .{ i, info.name });
        i += 1;
    }

    var buf: [10]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    try std.io.getStdIn().reader().streamUntilDelimiter(fbs.writer(), '\n', null);
    const input = std.mem.trim(u8, fbs.getWritten(), &std.ascii.whitespace);

    if (input.len > 0) {
        return try std.fmt.parseInt(c.PaDeviceIndex, input, 10) -| 1;
    } else if (default) |d| {
        return d;
    } else {
        return error.NoInputDevice;
    }
}

fn dump_bitstream() void {
    var bw = std.io.bufferedWriter(std.io.getStdOut().writer());

    while (true) {
        while (demodulated_data_read_pos != demodulated_data_write_pos) {
            const full = demodulated_data_write_pos + 1 == demodulated_data_read_pos or
                (demodulated_data_write_pos + 1 == demodulated_data_ring_buf.len and
                demodulated_data_read_pos == 0);
            if (full) {
                // For now, just drop this symbol and print a placeholder.
                bw.writer().writeByte('*') catch {};
            }

            @fence(.AcqRel); // TODO: use most appropriate ordering
            const byte = demodulated_data_ring_buf[demodulated_data_read_pos];
            bw.writer().print("{c}", .{byte}) catch {};

            @fence(.AcqRel); // TODO
            demodulated_data_read_pos = @intCast((demodulated_data_read_pos + 1) % demodulated_data_ring_buf.len);
        }

        bw.flush() catch {};
        std.time.sleep(1 * std.time.ns_per_ms);
    }
}

const Demodulator = modem.Demodulator(31, 44_100, 882);
var demodulator = Demodulator.init();

var demodulated_data_ring_buf: [1024]u8 = undefined;
var demodulated_data_write_pos: u16 = 0;
var demodulated_data_read_pos: u16 = 0;

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

    var sig = ins[0][0..frame_count];
    while (demodulator.demodulate(sig)) |res| {
        const read, const status = res;
        sig = sig[read..];
        switch (status) {
            .disconnected => produce_data('*'),
            .symbol => |sym| switch (sym) {
                .waiting, .ready => {},
                .payload => |p| produce_data(p),
            },
        }
    }

    return c.paContinue;
}

fn produce_data(byte: u8) void {
    @fence(.AcqRel); // TODO
    demodulated_data_ring_buf[demodulated_data_write_pos] = byte;

    @fence(.AcqRel); // TODO
    demodulated_data_write_pos = @intCast((demodulated_data_write_pos + 1) % demodulated_data_ring_buf.len);
}

fn paAssert(code: c.PaError) void {
    if (code != c.paNoError) {
        std.debug.panic("error: portaudio: {s}", .{c.Pa_GetErrorText(code)});
    }
}
