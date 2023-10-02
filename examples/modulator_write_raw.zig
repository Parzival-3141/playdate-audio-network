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

pub fn main() !void {
    const sample_rate: f32 = 44_100;
    const N = 30;

    var mod = modem.Modulator(N, sample_rate, 64).init();
    var demod = modem.Demodulator(N, sample_rate).init();

    const symbols = [_]?u8{
        'z',  'i',  'g',  '.',
        '.',  '.',  ' ',  null,
        null, null, null, null,
        null, null, null, null,
        'H',  'I',  ',',  ' ',
        'M',  'O',  'M',  '!',
        '\n', null, null, null,
        null, null, null, null,
    };

    const stdout = std.io.getStdOut();
    const stderr = std.io.getStdErr();
    for (symbols) |sym| {
        var buf: [N]f32 = undefined;
        mod.modulate(sym, &buf);
        for (buf) |float| {
            const bytes: [4]u8 = @bitCast(float);
            try stdout.writer().writeAll(&bytes);
        }

        const res = try demod.demodulate(&buf);
        if (res) |byte| {
            try stderr.writer().writeByte(byte);
            try stderr.writer().print("{c}\t0x{x:0<2}\n", .{ byte, byte });
        }
    }
}
