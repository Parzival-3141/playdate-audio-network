const std = @import("std");

/// Calculates the power of the given frequency within the input sample slice.
/// N is the number of DFT terms. Buffer must have at least N samples.
/// normalized_frequency is target frequency divided by sample rate.
/// Frequency detection is strongest when frequency is an integer multiple of (1/N * sample_rate).
///
/// Reference: https://en.wikipedia.org/wiki/Goertzel_algorithm#Power-spectrum_terms
pub fn goertzel(N: u16, buf: []const f32, normalized_frequency: f32) f32 {
    const coeff = 2 * @cos(2.0 * std.math.pi * normalized_frequency);

    var s_prev: f32 = 0;
    var s_prev2: f32 = 0;

    for (0..N) |i| {
        const samp = buf[i];
        const s = samp + (coeff * s_prev) - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }

    // Scaling the power by N squared makes it predictable for any value of N.
    const power = (s_prev2 * s_prev2) + (s_prev * s_prev) - (coeff * s_prev * s_prev2);
    const NN: f32 = @floatFromInt(N * N);
    return power / NN;
}

test goertzel {
    var buf: [100]f32 = undefined;
    try test_goertzel(&buf, 31, &.{.{ .k = 1 }}, 1, 0.25); // simple case
    try test_goertzel(&buf, 31, &.{.{ .k = 1, .phase = 0.321 }}, 1, 0.25); // initial phase has no effect
    try test_goertzel(&buf, 31, &.{.{ .k = 2 }}, 2, 0.25); // k has no effect
    try test_goertzel(&buf, 31, &.{.{ .k = 7 }}, 7, 0.25); // yet another k
    try test_goertzel(&buf, 70, &.{.{ .k = 1 }}, 1, 0.25); // N has no effect
    try test_goertzel(&buf, 70, &.{.{ .k = 33 }}, 33, 0.25); // different N and k
    try test_goertzel(&buf, 31, &.{.{ .k = 1, .amp = 0.1 }}, 1, 0.0025); // below unity gain
    try test_goertzel(&buf, 31, &.{.{ .k = 1, .amp = 4.0 }}, 1, 4.0); // above unity gain
    try test_goertzel(&buf, 100, &.{.{ .k = 33, .amp = 2.0 }}, 33, 1.0); // above unity gain, different N and k
    try test_goertzel(&buf, 31, &.{.{ .k = 5.1 }}, 5, 0.2366493); // signal is not perfectly k-aligned
    try test_goertzel(&buf, 31, &.{.{ .k = 5.5 }}, 5, 0.1069222); // signal is even less k-aligned
    try test_goertzel(&buf, 31, &.{.{ .k = 6.1 }}, 5, 0.001631582); // non-k-alinged non-target component causes spectral leakage
    try test_goertzel(&buf, 31, &.{.{ .k = 6 }}, 5, 0.0); // k-aligned non-target component goes by undetected
    try test_goertzel(&buf, 31, &.{ // several components present
        .{ .k = 1, .amp = 0.1, .phase = 0.1 },
        .{ .k = 2, .amp = 0.1, .phase = 0.2 },
        .{ .k = 3, .amp = 0.1, .phase = 0.3 },
        .{ .k = 4, .amp = 0.5, .phase = 0.4 }, // target
        .{ .k = 5, .amp = 0.1, .phase = 0.5 },
        .{ .k = 6, .amp = 0.1, .phase = 0.6 },
    }, 4, 0.0625);
    try test_goertzel(&buf, 31, &.{ // several non-k-aligned components present
        .{ .k = 0.987, .amp = 0.1, .phase = 0.1 },
        .{ .k = 2.1198, .amp = 0.1, .phase = 0.2 },
        .{ .k = 4, .amp = 0.5, .phase = 0.3 }, // target
        .{ .k = 6.42387, .amp = 0.1, .phase = 0.4 },
    }, 4, 0.0610713);
}

fn test_goertzel(
    buf: []f32,
    N: u16,
    sines: []const struct {
        /// k/N = normalized frequency for any sample rate
        k: f32,
        amp: f32 = 1,
        phase: f32 = 0,
    },
    /// Analyze this bin
    target_k: f32,
    mag: f32,
) !void {
    @memset(buf, 0);
    const N_float: f32 = @floatFromInt(N);
    for (sines) |sine| {
        var phase: f32 = sine.phase;
        for (buf) |*samp| {
            samp.* += @sin(phase * 2 * std.math.pi) * sine.amp;
            phase += sine.k / N_float;
        }
    }
    const target_freq = target_k / N_float;
    const res = goertzel(N, buf, target_freq);
    try std.testing.expectApproxEqAbs(mag, res, 0.00001);
}
