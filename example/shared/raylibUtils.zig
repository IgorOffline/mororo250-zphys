const math = @import("math");
const rl = @import("raylib");

/// Converts a math.Mat4x4 to a raylib Matrix
pub fn mathMat4ToRayLib(math_mat: math.Mat4x4) rl.Matrix {
    return rl.Matrix{
        .m0 = math_mat.v[0].v[0],
        .m1 = math_mat.v[0].v[1],
        .m2 = math_mat.v[0].v[2],
        .m3 = math_mat.v[0].v[3],
        .m4 = math_mat.v[1].v[0],
        .m5 = math_mat.v[1].v[1],
        .m6 = math_mat.v[1].v[2],
        .m7 = math_mat.v[1].v[3],
        .m8 = math_mat.v[2].v[0],
        .m9 = math_mat.v[2].v[1],
        .m10 = math_mat.v[2].v[2],
        .m11 = math_mat.v[2].v[3],
        .m12 = math_mat.v[3].v[0],
        .m13 = math_mat.v[3].v[1],
        .m14 = math_mat.v[3].v[2],
        .m15 = math_mat.v[3].v[3],
    };
}

/// Converts a raylib Color (u8 components 0-255) to a float array [4]f32 (0.0-1.0).
/// This conversion is necessary when passing color data to shaders via uniforms,
/// as GLSL shaders typically expect normalized float values.
pub fn colorToFloatArray(color: rl.Color) [4]f32 {
    return .{
        @as(f32, @floatFromInt(color.r)) / 255.0,
        @as(f32, @floatFromInt(color.g)) / 255.0,
        @as(f32, @floatFromInt(color.b)) / 255.0,
        @as(f32, @floatFromInt(color.a)) / 255.0,
    };
}
