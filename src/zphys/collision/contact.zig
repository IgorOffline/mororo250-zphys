const math = @import("math");

// todo: math.vec3 has a 16 bytes align improve
// Can we improve the memory footprint and performance of this struct somehow?
pub const Contact = struct {
    normal: math.Vec3,
    point_a: math.Vec3,
    point_b: math.Vec3,
    body_a: u32,
    body_b: u32,
    penetration: f32,
};

// todo decrease cache missed and improve memory footprint it might be a good idea to separate contact and points information
// Contact points a and contact points b are stored in world space
pub const ContactManifold = struct {
    body_a: u32,
    body_b: u33,
    normal: math.Vec3,
    penetration_depth: f32,
    contact_points_a: [4]math.Vec3,
    contact_points_b: [4]math.Vec3,
    length: u32,
};


/// Todo: Check paper constraints derivations for rigid body simulation in 3D - Daniel Chappuis -> U term
/// Jacobian constraint(-n, -(r_1 x n, n, r_2 x n)
/// n -> axis of the constraint
///
pub const PenetrationConstraint = struct {
    r1: math.Vec3, // lever arm 1
    r2: math.Vec3, // lever arm 2
    n: math.Vec3, // collision normal
    velocity_bias: math.Vec3, // Check slide 44 of Erin catto presentation: https://box2d.org/files/ErinCatto_ModelingAndSolvingConstraints_GDC2009.pdf -> treat bounce as velocity bias

    // Todo: should this really be stored or should we recalculate it?
    invert_inertia_n_x_r1: math.Vec3,
    invert_inertia_n_x_r2: math.Vec3,

    inverse_effective_mass: f32,
    accumulated_impulse: f32,
    inv_mass_a: f32,
    inv_mass_b: f32,
};
