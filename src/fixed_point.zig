const std = @import("std");

/// Numbers are scaled internally by 2^scale_bits.
pub fn Fixed(comptime Underlying: type, comptime scale_bits: u16) type {
    return struct {
        int: Underlying,

        const F = @This();

        pub fn from_float(x: f32) F {
            const factor: f32 = @floatFromInt(1 << scale_bits);
            const res: Underlying = @intFromFloat(@round(x * factor));
            return .{ .int = res };
        }

        pub fn to_float(x: F) f32 {
            const float: f32 = @floatFromInt(x.int);
            const factor: f32 = @floatFromInt(1 << scale_bits);
            return float / factor;
        }

        pub fn add(x: F, y: F) F {
            return .{ .int = x.int + y.int };
        }

        pub fn sub(x: F, y: F) F {
            return .{ .int = x.int - y.int };
        }

        pub fn mul(x: F, y: F) F {
            const res = (x.int * y.int) >> (scale_bits);
            return .{ .int = res };
        }

        // pub fn div(x: F, y: F) F {
        //     _ = y;
        //     _ = x;
        //     //
        // }
    };
}

const testing = std.testing;
test Fixed {
    const F = Fixed(u16, 14);
    const res = F.from_float(0.123).add(F.from_float(0.321));
    try testing.expectEqual(@as(u16, 7274), res.int);
    try testing.expectApproxEqAbs(@as(f32, 0.444), res.to_float(), 0.0001);
}
