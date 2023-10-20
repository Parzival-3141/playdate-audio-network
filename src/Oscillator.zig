//! Sine wave oscillator that uses a lookup table with linear interpolation.

const std = @import("std");
const assert = std.debug.assert;

const Oscillator = @This();

/// Current phase of the sine wave in range [0, 1).
phase: f32 = 0,

/// Produce the next oscillator amplitude and advance phase by phase_incr.
/// Phase increment can be calculated with (freq / sample_rate);
/// it is provided as a parameter instead of frequency because for this package because
/// it is more effecient to precalculate the phase_incrs of a finite set of frequencies.
///
/// phase_incr is asserted to be in the range [0, 0.5) to avoid crossing the Nyquist rate.
pub fn generate(osc: *Oscillator, sine_table: []const f32, phase_incr: f32) f32 {
    assert(phase_incr >= 0);
    assert(phase_incr < 0.5);
    assert(osc.phase >= 0 and osc.phase < 1);

    const out = sine_value(sine_table, osc.phase);

    osc.phase += phase_incr;
    while (osc.phase >= 1) osc.phase -= 1;

    return out;
}

/// Returns sin(phase * 2pi) using the lookup table with linear interpolation.
/// Phase must be in range [0, 1).
pub fn sine_value(sine_table: []const f32, phase: f32) f32 {
    assert(phase >= 0 and phase < 1);

    const index_float = phase * @as(f32, @floatFromInt(sine_table.len));
    const index: u32 = @intFromFloat(index_float);
    const frac = index_float - @as(f32, @floatFromInt(index));

    const y1 = sine_table[index % sine_table.len];
    const y2 = sine_table[(index + 1) % sine_table.len];

    return std.math.lerp(y1, y2, frac);
}

/// Creates a sinusoidal lookup table of the given length.
/// Holds amplitudes in range [-1, 1] for one cycle of a sine.
/// A length of 64 is fine for many practical purposes.
/// Lookup table is scaled ahead of time to compensate for expected amplitude summation at runtime.
pub fn create_sine_table(comptime len: u32) [len]f32 {
    if (len == 0) {
        @compileError("sine lookup table length must be > 0");
    }

    var vals: [len]f32 = undefined;
    for (&vals, 0..) |*v, i| {
        var frac: f64 = @floatFromInt(i);
        frac /= vals.len;
        v.* = @floatCast(@sin(frac * std.math.pi * 2));
    }

    return vals;
}
