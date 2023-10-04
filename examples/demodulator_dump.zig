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

const Demodulator = @import("modem").modem2.Demodulator;

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
        .suggestedLatency = in_info.*.defaultHighInputLatency,
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

const goertzel_N = 96;
const resync_adjust = 0;

const Demod = Demodulator(goertzel_N, sample_rate);
var demodulator = Demod.init();

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
    _ = frame_count;
    _ = out_buf;
    _ = user_data;
    _ = status_flags;
    _ = time_info;

    const ins: [*]const [*]const f32 = if (in_buf) |buf| @alignCast(@ptrCast(buf)) else unreachable;

    var i: u32 = 0;
    _ = i;

    var offset: u32 = 0;
    while (demodulator.demodulate(ins[0][offset..])) |res| : (offset += Demod.window_len) switch (res) {
        .no_signal => produce_data('*'),
        .not_ready, .ready => {},
        .payload => |byte| produce_data(byte),
    };

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
