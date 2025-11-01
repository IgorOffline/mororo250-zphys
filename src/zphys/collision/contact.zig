const math = @import("math");

// todo: math.vec3 has a 16 bytes align improve
// Can we improve the memory footprint and performance of this struct somehow?
pub const Contact = struct {
    normal: math.Vec3,
    point_local_a: math.Vec3,
    point_local_b: math.Vec3,
    body_a: u32,
    body_b: u32,
    penetration: f32,
};

// todo decrease cache missed and improve memory footprint it might be a good idea to separate contact and points information
pub const ContactManifold = struct {
    body_a: u32,
    body_b: u32,
    normal: math.vec3,
    penetration_depth: f32,
    contact_points_a: [4]math.vec3,
    contact_points_b: [4]math.vec3,
    length: u32,
};
