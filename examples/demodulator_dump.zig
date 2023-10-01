const std = @import("std");
const dbprint = std.debug.print;
const c = @cImport({
    @cInclude("portaudio.h");
});

const demodulator = @import("modem").demodulator;

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
        44_100,
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

const goertzel_N = 36;
/// Duration in samples of a single symbol. Inversely proportional to baud.
const goertzel_symbol_period = 36;
var sample_counter: u32 = 0;
var goertzel_buffer: [goertzel_symbol_period]f32 = undefined;

const demodulator_base_freq = sample_rate / goertzel_N;
const Demod = demodulator.Demodulator(1, &.{
    .{ demodulator_base_freq * 1, demodulator_base_freq * 2 },
    .{ demodulator_base_freq * 3, demodulator_base_freq * 4 },
    .{ demodulator_base_freq * 5, demodulator_base_freq * 6 },
    .{ demodulator_base_freq * 7, demodulator_base_freq * 8 },
    .{ demodulator_base_freq * 9, demodulator_base_freq * 10 },
    .{ demodulator_base_freq * 11, demodulator_base_freq * 12 },
    .{ demodulator_base_freq * 13, demodulator_base_freq * 14 },
    .{ demodulator_base_freq * 15, demodulator_base_freq * 16 },
}, sample_rate, goertzel_N);
var demod = Demod{};

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

    for (0..frame_count) |i| {
        demodulate_stream(ins[0][i]);
    }

    return c.paContinue;
}

fn demodulate_stream(in_samp: f32) void {
    goertzel_buffer[sample_counter] = in_samp;
    sample_counter += 1;

    if (sample_counter == goertzel_buffer.len) {
        sample_counter = 0;

        @fence(.AcqRel); // TODO
        const symbol = demod.analyze(&goertzel_buffer);
        demodulated_data_ring_buf[demodulated_data_write_pos] = @bitCast(symbol);

        @fence(.AcqRel); // TODO
        demodulated_data_write_pos = @intCast((demodulated_data_write_pos + 1) % demodulated_data_ring_buf.len);
    }
}

fn paAssert(code: c.PaError) void {
    if (code != c.paNoError) {
        std.debug.panic("error: portaudio: {s}", .{c.Pa_GetErrorText(code)});
    }
}
