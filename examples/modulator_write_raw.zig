const std = @import("std");
var rng = std.rand.DefaultPrng.init(12345678);

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
    const N = 35;

    var mod = modem.Modulator(N, sample_rate, 64).init();
    var demod = modem.Demodulator(N, sample_rate).init();

    const Message = union(enum) {
        nulls: u16,
        str: []const u8,
    };
    const messages = [_]Message{
        .{ .str = "zig... " },
        .{ .nulls = 7 },
        .{ .str = "HI MOM!\n" },
        .{ .nulls = 15 },
        .{ .str = "ok; " },
        .{ .str = "Insert your favorite song or poem here: " },
        .{ .nulls = 1 },
        .{ .str = "0123456789 ABCDEFGHIJKLMNOPQRSTUVWXYZ\n" },
        .{ .nulls = 3 },
        .{ .str = 
        \\One morning, as Gregor Samsa was waking up from anxious dreams, he
        \\discovered that in bed he had been changed into a monstrous vermin.
        \\He lay on his armour-hard back and saw, as he lifted his head up a little,
        \\his brown, arched abdomen divided up into rigid bow-like sections.
        \\From this height the blanket, just about ready to slide off completely,
        \\could hardly stay in place. His numerous legs, pitifully thin in comparis-
        \\on to the rest of his circumference, flickered helplessly before his eyes.
        \\
        \\
        },
    };
    // const messages = [1]Message{.{ .str = &.{ 0, 1, 2, 3, 4, 5, 6, 7 } }} ** 100;

    const symbols = comptime blk: {
        comptime var syms: []const ?u8 = &.{};
        inline for (messages) |msg| {
            switch (msg) {
                .str => |str| {
                    for (str) |char| syms = syms ++ &[_]?u8{char};
                },
                .nulls => |count| {
                    for (0..count) |_| syms = syms ++ &[_]?u8{null};
                },
            }
        }
        break :blk syms;
    };

    const stdout = std.io.getStdOut();
    const stderr = std.io.getStdErr();

    var buf: [N * 2]f32 = undefined;
    var misses: u32 = 0;

    // Change this to simulate more/less sync drift.
    // const sync_error = 0.15;
    const sync_error = 0.04;
    const offset: u32 = @intFromFloat(sync_error * N);

    // Change this to add more or less white noise to the signal.
    const noise_floor = 0.01;
    // const noise_floor = 0.025;

    std.log.info("applying sync error: {d:.4} ({} samples)", .{ sync_error, offset });
    std.log.info("applying noise floor: {d:.4}", .{noise_floor});
    std.log.info("relaying signal...", .{});

    for (symbols, 0..) |sym, i| {
        // Make a buffer that can hold two periods,
        // and simulate out of sync symbol windows.
        if (i == 0) {
            mod.modulate(sym, buf[0..N]);
            mod.modulate(sym, buf[N..]);
        } else {
            std.mem.copyBackwards(f32, buf[0..N], buf[N..]);
            mod.modulate(sym, buf[N..]);
        }

        // Add some noise to the signal
        for (&buf) |*samp| {
            const noise_val = rng.random().float(f32) * noise_floor * 2.0; // [0, 2*noise]
            samp.* += noise_val - (noise_val / 2.0); // [-noise/2, noise/2]
        }

        const offset_buf = buf[N - offset ..][0..N];

        for (offset_buf) |float| {
            const bytes: [4]u8 = @bitCast(float);
            try stdout.writer().writeAll(&bytes);
        }

        // const res = try demod.demodulate(offset_buf);
        const res = demod.demodulate(offset_buf) catch blk: {
            try stderr.writer().writeByte('*');
            misses += 1;
            break :blk null;
        };
        if (res) |byte| {
            try stderr.writer().writeByte(byte);
            // try stderr.writer().print("{c}\t0x{x:0>2}\n", .{ byte, byte });
        }
    }

    if (misses > 0) {
        std.log.err("total misses: {}", .{misses});
        std.os.exit(1);
    } else {
        std.log.info("no misses", .{});
    }
}
