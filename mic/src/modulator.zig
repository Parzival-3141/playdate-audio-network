const std = @import("std");

/// Create a multiple frequency-shift keying (MFSK) modulator.
///
/// The modulator stores `osc_freqs` number of oscillators, each of which
/// selects one from a set of frequencies to encode the currently stored symbol
/// as a single sine wave; the total number of sine waves produced is `osc_freqs.len`.
/// `osc_freqs` specifies the selectable frequencies for each oscillator.
///
/// For practical purposes, all frequencies should be unique and disparate enough to be
/// clearly detectable with signal analysis. Further, they should ideally be at the center
/// of their respective DFT frequency bin.
///
/// Because MFSK requires 2^(bits_per_oscillator) unique frequencies for each oscillator,
/// small values of bits_per_oscillator may be preferred.
///
/// The current symbol being encoded is represented as a packed struct.
/// Each field represents some of the bits of the data, and simultaneously the selected
/// frequency index for the corresponding oscillator.
///
/// The oscillators use a sinusoidal lookup table with linear interpolation.
/// 64 is a reasonable table size for practical purposes.
pub fn Modulator(
    comptime bits_per_oscillator: u4,
    comptime osc_freqs: []const [1 << bits_per_oscillator]f32,
    comptime sample_rate: f32,
    comptime sine_table_len: u32,
) type {
    if (osc_freqs.len == 0) {
        @compileError("Modulator requires at least one oscillator");
    }

    const freqs_per_osc = 1 << bits_per_oscillator;

    // Convert frequencies into phase increments. Since we have a fixed set of
    // selectable frequencies, precalculating these should be more efficient than
    // doing floating point division every time we generate a sample value.
    comptime var phase_incrs: [osc_freqs.len * freqs_per_osc]f32 = undefined;
    inline for (osc_freqs, 0..) |osc, i| {
        for (osc, 0..) |freq, j| {
            phase_incrs[i * freqs_per_osc + j] = freq / sample_rate;
        }
    }

    // Construct a packed struct where each field holds the bits each oscillator.
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
        const Mod = @This();

        pub const SymbolInt = @Type(.{
            .Int = .{
                .signedness = .unsigned,
                .bits = osc_freqs.len * bits_per_oscillator,
            },
        });

        /// Symbol is a packed struct containing all the bits encoded by the
        /// selected or detected freqeuncies within one messaging period.
        /// It is a packed struct where each field represents the current frequency index
        /// of a single oscillator.
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

        pub const sine_table = create_lookup_table(sine_table_len);

        oscillators: [osc_freqs.len]Oscillator,
        /// Current symbol value. Each field's bits also equals the corresponding
        /// oscillator's frequency index.
        symbol: Symbol,

        pub fn init() Mod {
            var mod = Mod{
                .oscillators = [1]Oscillator{.{}} ** osc_freqs.len,
                .symbol = undefined,
            };
            inline for (@typeInfo(Symbol).Struct.fields) |field| {
                @field(mod.symbol, field.name) = 0;
            }
            return mod;
        }

        /// For each oscillator, produce the next amplitude value and update phase.
        pub fn generate_sample(mod: *Mod) [osc_freqs.len]f32 {
            var out: [osc_freqs.len]f32 = undefined;
            inline for (&mod.oscillators, @typeInfo(Symbol).Struct.fields, 0..) |*osc, field, i| {
                const val = @field(mod.symbol, field.name);
                out[i] = osc.generate(&sine_table, phase_incrs[i * freqs_per_osc + val]);
            }
            return out;
        }
    };
}

/// Creates a sinusoidal lookup table of the given length.
/// Holds amplitudes in range [-1, 1] for one cycle of a sine.
/// A length of 64 is fine for many practical purposes.
pub fn create_lookup_table(comptime len: u32) [len]f32 {
    if (len == 0) {
        @compileError("sine lookup table length must be > 0");
    }

    var vals: [len]f32 = undefined;
    for (&vals, 0..) |*v, i| {
        v.* = @floatCast(@sin(i / @as(f64, vals.len) * std.math.pi * 2));
    }

    return vals;
}

/// Sine wave oscillator that uses a lookup table with linear interpolation.
pub const Oscillator = struct {
    /// Current phase of the sine wave in range [0, 1).
    phase: f32 = 0,

    /// Produce the next oscillator amplitude and advance phase by phase_incr.
    /// Phase increment can be calculated with (freq / sample_rate);
    /// it is provided as a parameter instead of frequency because for this package because
    /// it is more effecient to precalculate the phase_incrs of a finite set of frequencies.
    ///
    /// phase_incr is asserted to be positive.
    pub fn generate(osc: *Oscillator, comptime sine_table: []const f32, phase_incr: f32) f32 {
        std.debug.assert(phase_incr >= 0);

        const index_float = osc.phase * sine_table.len;
        const index: u32 = @intFromFloat(index_float);
        const frac = index_float - @as(f32, @floatFromInt(index));

        const y1 = sine_table[index % sine_table.len];
        const y2 = sine_table[(index + 1) % sine_table.len];
        const out = y1 * (1 - frac) + (y2 * frac);

        osc.phase += phase_incr;
        while (osc.phase >= 1) osc.phase -= 1;

        return out;
    }
};

test Modulator {
    const Mod = Modulator(2, &.{
        .{ 100, 110, 120, 130 },
        .{ 200, 210, 220, 230 },
        .{ 300, 310, 320, 330 },
        .{ 400, 410, 420, 430 },
    }, 44_100, 64);
    var mod = Mod.init();

    try std.testing.expect(Mod.SymbolInt == u8);
    try std.testing.expect(@bitSizeOf(Mod.Symbol) == 8);
    try std.testing.expect(@TypeOf(mod.symbol.osc0) == u2);
    try std.testing.expect(@TypeOf(mod.symbol.osc1) == u2);
    try std.testing.expect(@TypeOf(mod.symbol.osc2) == u2);
    try std.testing.expect(@TypeOf(mod.symbol.osc3) == u2);

    for ("Hello, world!") |byte| {
        mod.symbol = @bitCast(byte);
        _ = mod.generate_sample();
    }
}
