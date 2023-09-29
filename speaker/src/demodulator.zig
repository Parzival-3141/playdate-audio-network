const std = @import("std");

pub fn Demodulator(
    comptime bits_per_oscillator: u4,
    comptime osc_freqs: []const [1 << bits_per_oscillator]f32,
    comptime sample_rate: f32,
    comptime goertzel_N: u16,
) type {
    const freqs_per_osc = 1 << bits_per_oscillator;

    comptime var normalized_freqs: [osc_freqs.len][freqs_per_osc]f32 = undefined;
    for (osc_freqs, 0..) |osc, i| {
        for (osc, 0..) |freq, j| {
            normalized_freqs[i][j] = freq / sample_rate;
        }
    }

    comptime var symbol_fields: [osc_freqs.len]std.builtin.Type.StructField = undefined;
    inline for (0..osc_freqs.len) |i| {
        const Int = @Type(.{
            .Int = .{
                .signedness = .unsigned,
                .bits = bits_per_oscillator,
            },
        });

        var buf: [10]u8 = undefined;
        const i_str = comptime try std.fmt.bufPrint(&buf, "osc{}", .{i});

        symbol_fields[i] = .{
            .name = i_str,
            .type = Int,
            .default_value = null,
            .is_comptime = false,
            .alignment = 0,
        };
    }

    return struct {
        const Demod = @This();

        pub const SymbolInt = @Type(.{
            .Int = .{
                .signedness = .unsigned,
                .bits = osc_freqs.len * bits_per_oscillator,
            },
        });

        pub const Symbol = @Type(.{
            .Struct = .{
                .layout = .Packed,
                .backing_integer = SymbolInt,
                .fields = &symbol_fields,
                .decls = &[0]std.builtin.Type.Declaration{},
                // TODO: Why were packed tuples removed?
                // https://github.com/ziglang/zig/pull/16551
                .is_tuple = false,
            },
        });

        pub fn analyze(demod: Demod, buf: []f32) Symbol {
            _ = demod;

            var symbol: Symbol = undefined;

            const Int = @Type(.{
                .Int = .{
                    .signedness = .unsigned,
                    .bits = bits_per_oscillator,
                },
            });

            inline for (normalized_freqs, 0..) |osc, i| {
                var strongest_index: Int = 0;
                var strongest_mag: f32 = 0;
                for (osc, 0..) |freq, j| {
                    const mag = goertzel(buf, freq, goertzel_N);
                    if (mag > strongest_mag) {
                        strongest_mag = mag;
                        strongest_index = @intCast(j);
                    }
                }

                @field(symbol, @typeInfo(Symbol).Struct.fields[i].name) = strongest_index;
            }

            return symbol;
        }
    };
}

/// Calculates the power of the given frequency within the input sample slice.
/// N is the number of DFT terms.
/// normalized_frequency is target frequency divided by sample rate.
/// Frequency detection is strongest when frequency is an integer multiple of (1/N * sample_rate).
///
/// Reference: https://en.wikipedia.org/wiki/Goertzel_algorithm#Power-spectrum_terms
pub fn goertzel(buf: []const f32, normalized_frequency: f32, comptime N: u16) f32 {
    const coeff = 2 * @cos(2.0 * std.math.pi * normalized_frequency);

    var s_prev: f32 = 0;
    var s_prev2: f32 = 0;

    for (0..N) |i| {
        const s = buf[i] + (coeff * s_prev) - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }

    return (s_prev2 * s_prev2) + (s_prev * s_prev) - (coeff * s_prev * s_prev2);
}
