const std = @import("std");

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

const sample_rate: f32 = 44_100;
const N = 31;
const baud = 882;

const Modulator = modem.Modulator(N, sample_rate, baud, 64);
const Demodulator = modem.Demodulator(N, sample_rate, baud);
var mod = Modulator.init();
var demod = Demodulator.init();

var sample_buf: [Modulator.symbol_len]f32 = undefined;

pub fn main() !void {
    const symbols = [_]modem.Symbol{
        .waiting,            .waiting,            .waiting,            .ready,
        .{ .payload = 't' }, .{ .payload = 'a' }, .{ .payload = 'c' }, .{ .payload = 'o' },
        .{ .payload = 'h' }, .{ .payload = 'a' }, .{ .payload = 't' }, .ready,
        .ready,              .ready,              .ready,              .ready,
    };

    for (symbols) |symbol| {
        mod.modulate(symbol, &sample_buf);
        try write_out_samples(&sample_buf);
    }

    std.debug.print("\n", .{});
    if (disconnects > 0) {
        std.log.err("got {d} disconnected windows", .{disconnects});
        std.os.exit(1);
    }
}

var disconnects: usize = 0;

fn write_out_samples(signal: []const f32) !void {
    const stderr = std.io.getStdErr().writer();
    const stdout = std.io.getStdOut().writer();

    for (signal) |float| {
        const buf_bytes: [4]u8 = @bitCast(float);
        try stdout.writeAll(&buf_bytes);
    }

    var sig = signal;
    while (demod.demodulate(signal)) |res| {
        const read, const status = res;
        sig = sig[read..];
        switch (status) {
            .disconnected => {
                try stderr.writeByte('*');
                disconnects += 1;
            },
            .symbol => |symbol| switch (symbol) {
                .waiting, .ready => {},
                .payload => |byte| try stderr.writeByte(byte),
            },
        }
    }
}
