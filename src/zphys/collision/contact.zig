const math = @import("math");

// todo: math.vec3 has a 16 bytes align improve
// Can we improve the memory footprint and performance of this struct somehow?
pub const Contact = struct {
    normal: math.Vec3,
    contact_point_a: math.Vec3,
    contact_point_b: math.Vec3,
    body_a: u32,
    body_b: u32,
    penetration: f32,
    friction: f32,
    restitution: f32,
};

pub const ContactManifold = struct {
    body_a: u32,
    boody_b: u32,
    normal: math.vec3,
    penetration_depth: f32,

    // to decrease cache missed and improve memory footprint it might be a good idea to separate contact and points information
    relative_contact_poinit_on_a: []math.vec3,
    relative_contact_point_on_b: []math.vec3,
};
