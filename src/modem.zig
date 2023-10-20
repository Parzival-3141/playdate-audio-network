pub const Modulator = @import("modulator.zig").Modulator;
pub const Demodulator = @import("demodulator.zig").Demodulator;
pub const goertzel = @import("goertzel.zig").goertzel;
pub const Oscillator = @import("Oscillator.zig");

/// Values that one side (client) can send to the other (server).
/// Both sides may simultaneously act as client and server to one another.
pub const Symbol = union(enum) {
    /// Client is running but has not detected a valid signal from the server.
    /// Server should not send a payload back while the client is in this state.
    waiting,
    /// Client has detected the server but there is no data to send.
    ready,
    /// Client has detected the server and is sending a byte of data.
    payload: u8,
};

pub const OscillatorRates = extern struct {
    header: [3]f32,
    payload: [8]f32,

    pub const osc_count = 11;

    /// Returns a set of frequencies, normalized to a fraction [0, 1] of sample rate,
    /// that are optimal for Goertzel detection with N DFT terms.
    ///
    /// One symbol requires 11 unique frequencies, in increasing order:
    /// - Three for the header oscillator, waiting/ready/payload (FSK)
    /// - Eight for each payload bit oscillator, in least to most signifcant order (ASK)
    pub fn init(N: u16, sample_rate: f32) OscillatorRates {
        const base_freq = sample_rate / @as(f32, @floatFromInt(N));
        if (base_freq < 20) {
            @compileError("modem signal contains a frequency below 20 Hz; choose a smaller value for N");
        }
        if (base_freq * osc_count > 20_000) {
            @compileError("modem signal contains a frequency over 20 kHz; choose a larger value for N");
        }

        var rates: [osc_count]f32 = undefined;
        for (&rates, 0..) |*f, i| {
            f.* = (base_freq * (i + 1)) / sample_rate;
        }

        return @bitCast(rates);
    }
};

test "modulate and demodulate symbols" {
    try test_modulate_and_demodulate(26, 1_000, 10);
    try test_modulate_and_demodulate(26, 44_100, 100);
    try test_modulate_and_demodulate(26, 44_100, 441);
    try test_modulate_and_demodulate(30, 48_000, 500);
    try test_modulate_and_demodulate(40, 10_000, 50);
    try test_modulate_and_demodulate(40, 14_000, 160);
    try test_modulate_and_demodulate(60, 48_000, 10);
    try test_modulate_and_demodulate(60, 48_000, 50);
}

const std = @import("std");

fn test_modulate_and_demodulate(
    comptime N: u16,
    comptime sample_rate: comptime_float,
    comptime baud: comptime_float,
) !void {
    const Mod = Modulator(N, sample_rate, baud, 64);
    const Demod = Demodulator(N, sample_rate, baud);

    var mod = Mod.init();
    var demod = Demod.init();

    const symbols = [_]Symbol{
        .waiting,            .waiting,            .waiting,            .ready,
        .{ .payload = 't' }, .{ .payload = 'a' }, .{ .payload = 'c' }, .{ .payload = 'o' },
        .{ .payload = 'h' }, .{ .payload = 'a' }, .{ .payload = 't' }, .ready,
        .ready,              .ready,              .ready,              .ready,
    };

    var signal: [Mod.symbol_len]f32 = undefined;
    var results: [symbols.len]Symbol = undefined;
    var results_index: usize = 0;

    for (symbols, 0..) |sym, i| {
        mod.modulate(sym, &signal);

        var pos: usize = 0;
        while (pos < signal.len) {
            const read, const status = demod.demodulate(signal[pos..]) orelse break;
            pos += read;
            switch (status) {
                .disconnected => {
                    std.debug.panic("could not demodulate; iter {d}, expected {}", .{ i, sym });
                },
                .symbol => |s| {
                    results[results_index] = s;
                    results_index += 1;
                },
            }
        }
    }

    try std.testing.expectEqual(symbols.len, results_index);
    try std.testing.expectEqualSlices(Symbol, &symbols, &results);
}
