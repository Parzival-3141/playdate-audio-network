const std = @import("std");
const assert = std.debug.assert;

const modem = @import("modem.zig");
const OscillatorRates = modem.OscillatorRates;
const Symbol = modem.Symbol;
const Oscillator = @import("Oscillator.zig");

pub const Options = struct {
    sine_table_len: u32 = 64,
};

/// Creates a type that modulates data as an audio signal.
/// The signal uses a combination of frequency-shift keying (FSK) and amplitude-shifk keying (ASK).
/// Sine waves are generated using a lookup table of the given length.
/// 64 is a reasonable table length.
pub fn Modulator(
    comptime N: u16,
    comptime sample_rate: comptime_float,
    comptime baud: comptime_float,
    comptime opts: Options,
) type {
    comptime {
        assert(N > 0);
        assert(sample_rate > 0);
        assert(baud > 0);
    }

    const rates = OscillatorRates.init(N, sample_rate);

    return struct {
        const Mod = @This();

        const symbol_len_float = @floor(sample_rate / baud);
        pub const symbol_len: u16 = @intFromFloat(symbol_len_float);
        comptime {
            assert(symbol_len == symbol_len_float);
            if (N > symbol_len) {
                @compileError("symbol period must be greater than N samples");
            }
        }

        pub const sine_table = Oscillator.create_sine_table(opts.sine_table_len);

        header_osc: Oscillator,
        payload_oscs: [8]Oscillator,

        pub fn init() Mod {
            return Mod{
                .header_osc = .{ .phase = 0 },
                .payload_oscs = .{
                    .{ .phase = 0.1 },
                    .{ .phase = 0.2 },
                    .{ .phase = 0.3 },
                    .{ .phase = 0.4 },
                    .{ .phase = 0.5 },
                    .{ .phase = 0.6 },
                    .{ .phase = 0.7 },
                    .{ .phase = 0.8 },
                },
            };
        }

        /// Generate signal for the symbol into buf.
        /// Buf must be large enough to store one symbol period.
        pub fn modulate(mod: *Mod, symbol: Symbol, buf: []f32) void {
            const header = @intFromEnum(symbol);
            for (buf[0..symbol_len]) |*out| {
                var samp = mod.header_osc.generate(&sine_table, rates.header[header]);
                out.* = samp;
            }

            const byte = switch (symbol) {
                .payload => |b| b,
                else => 0,
            };

            for (0..8) |bit| {
                const mask: u8 = @as(u8, 1) << @intCast(bit);
                const bit_set = byte & mask != 0;
                if (!bit_set) {
                    continue;
                }

                for (buf[0..symbol_len]) |*out| {
                    var samp = mod.payload_oscs[bit].generate(&sine_table, rates.payload[bit]);
                    out.* += samp;
                }
            }

            // Total scale factor is inverse of max number of
            // simultaneously mixed oscillators to avoid clipping
            const scale: f32 = 1.0 / @as(f32, OscillatorRates.osc_count);

            // Shape symbol window with a full sine cycle.
            for (buf[0..symbol_len], 0..) |*out, i| {
                var frac = @as(f32, @floatFromInt(i)) / symbol_len_float + 0.75;
                if (frac >= 1) frac -= 1;
                var amp = Oscillator.sine_value(&sine_table, frac); // [-1, 1]
                amp *= 0.5; // [-0.5, 0.5]
                amp += 0.5; // [0, 1]
                amp *= scale;
                out.* *= amp;
            }
        }
    };
}

test Modulator {
    const M = Modulator(26, 44100, 441, .{});
    var m = M.init();
    var buf: [100]f32 = undefined;
    m.modulate(.waiting, &buf);
    m.modulate(.waiting, &buf);
    m.modulate(.ready, &buf);
    m.modulate(.{ .payload = 'h' }, &buf);
    m.modulate(.{ .payload = 'i' }, &buf);
    m.modulate(.ready, &buf);
}
