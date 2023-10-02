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

const goertzel_N = 30;
const resync_adjust = -35;

var demodulator = Demodulator(goertzel_N, sample_rate).init();

var demodulated_data_ring_buf: [1024]u8 = undefined;
var demodulated_data_write_pos: u16 = 0;
var demodulated_data_read_pos: u16 = 0;

var demodulate_signal_buf_arr: [goertzel_N]f32 = undefined;
var demodulate_signal_buf: []f32 = demodulate_signal_buf_arr[0..0];

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

    var i: u32 = 0;

    // Left over from the previous callback
    if (demodulate_signal_buf.len > 0) {
        if (frame_count < goertzel_N) @panic("TODO: handle small callback windows");

        while (i + goertzel_N < demodulate_signal_buf.len) : (i += goertzel_N) {
            const res = demodulator.demodulate(demodulate_signal_buf[i .. i + goertzel_N]) catch no_sig: {
                // Add delay to try and resync
                const adjusted = @as(i64, @intCast(i)) + resync_adjust;
                if (adjusted > 0 and adjusted < demodulate_signal_buf.len) {
                    i = @intCast(adjusted);
                }
                break :no_sig null;
            };
            if (res) |char| produce_data(char);
        }

        demodulate_signal_buf = demodulate_signal_buf_arr[0..0];
    }

    // Incoming stream data
    while (i + goertzel_N < frame_count) : (i += goertzel_N) {
        const res = demodulator.demodulate(ins[0][i .. i + goertzel_N]) catch no_sig: {
            // Add delay to try and resync
            const adjusted = @as(i64, @intCast(i)) + resync_adjust;
            if (adjusted > 0 and adjusted < frame_count) {
                i = @intCast(adjusted);
            }
            break :no_sig null;
        };
        if (res) |char| produce_data(char);
    }

    // There are left over samples to buffer for next callback
    if (i < frame_count) {
        const rem = frame_count - i;
        @memcpy(demodulate_signal_buf_arr[0..rem], ins[0][i..frame_count]);
        demodulate_signal_buf = demodulate_signal_buf_arr[0..rem];
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
