const std = @import("std");
const math = @import("math");

/// Compute barycentric coordinates of point `p` with respect to the line segment AB.
/// Returns weights (u, v) such that p ≈ u*A + v*B and u + v = 1.
/// - Degenerate case (A == B): returns (1, 0).
pub inline fn barycentricLine(p: math.Vec3, a: math.Vec3, b: math.Vec3) math.Vec2 {
    const ab = b.sub(&a);
    const denom = ab.dot(&ab);
    if (denom <= math.eps_f32) {
        // Degenerate line: fully weight A
        return math.vec2(1.0, 0.0);
    }
    const t = p.sub(&a).dot(&ab) / denom;
    return math.vec2(1.0 - t, t);
}

// Todo: Optimize to use matrix multiplication instead
/// Compute barycentric coordinates of point `p` with respect to triangle ABC.
/// Returns weights (u, v, w) such that p ≈ u*A + v*B + w*C and u + v + w = 1.
/// - Degenerate case (area ~ 0): falls back to the longest edge and uses line barycentric,
///   mapping the resulting 2D weights into 3D (the opposite vertex gets weight 0).
pub inline fn barycentricTriangle(p: math.Vec3, a: math.Vec3, b: math.Vec3, c: math.Vec3) math.Vec3 {
    const v0 = b.sub(&a);
    const v1 = c.sub(&a);
    const v2 = p.sub(&a);

    const d00 = v0.dot(&v0);
    const d01 = v0.dot(&v1);
    const d11 = v1.dot(&v1);
    const d20 = v2.dot(&v0);
    const d21 = v2.dot(&v1);

    const denom = d00 * d11 - d01 * d01;

    // If triangle is degenerate (area ~ 0), use the longest edge fallback.
    if (denom <= math.eps_f32) {
        const ab2 = d00;
        const ac2 = d11;
        const bc = c.sub(&b);
        const bc2 = bc.dot(&bc);

        if (ab2 >= ac2 and ab2 >= bc2) {
            const w2 = barycentricLine(p, a, b);
            return math.vec3(w2.x(), w2.y(), 0.0);
        } else if (ac2 >= ab2 and ac2 >= bc2) {
            const w2 = barycentricLine(p, a, c);
            return math.vec3(w2.x(), 0.0, w2.y());
        } else {
            const w2 = barycentricLine(p, b, c);
            return math.vec3(0.0, w2.x(), w2.y());
        }
    }

    const v = (d11 * d20 - d01 * d21) / denom;
    const w = (d00 * d21 - d01 * d20) / denom;
    const u = 1.0 - v - w;
    return math.vec3(u, v, w);
}

test "barycentricLine basic" {
    const a = math.vec3(0.0, 0.0, 0.0);
    const b = math.vec3(1.0, 0.0, 0.0);
    const p = math.vec3(0.25, 0.0, 0.0);

    const w = barycentricLine(p, a, b);

    try std.testing.expect(math.eql(f32, w.x, 0.75, 1e-6));
    try std.testing.expect(math.eql(f32, w.y, 0.25, 1e-6));
    // Sum to 1
    try std.testing.expect(math.eql(f32, w.x + w.y, 1.0, 1e-6));
}

test "barycentricTriangle basic" {
    const a = math.vec3(0.0, 0.0, 0.0);
    const b = math.vec3(1.0, 0.0, 0.0);
    const c = math.vec3(0.0, 1.0, 0.0);
    const p = math.vec3(0.25, 0.25, 0.0);

    const w = barycentricTriangle(p, a, b, c);

    try std.testing.expect(math.eql(f32, w.x, 0.5, 1e-6));
    try std.testing.expect(math.eql(f32, w.y, 0.25, 1e-6));
    try std.testing.expect(math.eql(f32, w.z, 0.25, 1e-6));
    // Sum to 1
    try std.testing.expect(math.eql(f32, w.x + w.y + w.z, 1.0, 1e-6));
}

test "barycentricTriangle degenerate fallback to longest edge" {
    const a = math.vec3(0.0, 0.0, 0.0);
    const b = math.vec3(1.0, 0.0, 0.0);
    const c = math.vec3(2.0, 0.0, 0.0); // Colinear (degenerate triangle), longest edge is AC
    const p = math.vec3(0.5, 0.0, 0.0);

    const w = barycentricTriangle(p, a, b, c);

    // Longest edge AC => line weights (0.75, 0.25) mapped to (u, v, w) = (0.75, 0.0, 0.25)
    try std.testing.expect(math.eql(f32, w.x, 0.75, 1e-6));
    try std.testing.expect(math.eql(f32, w.y, 0.0, 1e-6));
    try std.testing.expect(math.eql(f32, w.z, 0.25, 1e-6));
    try std.testing.expect(math.eql(f32, w.x + w.y + w.z, 1.0, 1e-6));
}
