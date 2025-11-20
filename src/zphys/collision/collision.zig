const std = @import("std");
const math = @import("math");
const body_module = @import("../body.zig");
const MotionComp = body_module.MotionComp;
const TransformComp = body_module.TransformComp;
const PhysicsPropsComp = body_module.PhysicsPropsComp;
const BodyComponents = body_module.BodyComponents;
const Shape = @import("shape.zig").Shape;

const contact = @import("contact.zig");
const sphere_sphere = @import("sphere_sphere.zig");
const sphere_box = @import("sphere_box.zig");
const box_box = @import("box_box.zig");
const gjk = @import("gjk.zig");
const sat = @import("sat.zig");

pub const Contact = contact.Contact;

pub const collideSphereSphere = sphere_sphere.collideSphereSphere;
pub const collideSphereBox = sphere_box.collideSphereBox;
pub const collideBoxBox = box_box.collideBoxBox;

pub const satBoxBoxContact = sat.satBoxBoxContact;
pub const ContactManifold = contact.ContactManifold;

// Todo: Add BroadPhase collision check in here
pub fn generateContacts(
    bodies: std.MultiArrayList(BodyComponents).Slice,
    contacts_out: *std.ArrayList(Contact),
    manifolds_out: *std.ArrayList(contact.ContactManifold)
) void {
    const transforms = bodies.items(.transform);
    const shapes = bodies.items(.shape);
    const physics_props = bodies.items(.physics_props);
    
    var index_a: usize = 0;
    while (index_a < bodies.len) : (index_a += 1) {
        var index_b: usize = index_a + 1;
        while (index_b < bodies.len) : (index_b += 1) {
            // Skip if both bodies are static
            // After separating statics
            if (physics_props[index_a].inverseMass == 0 and physics_props[index_b].inverseMass == 0) continue;

            switch (shapes[index_a]) {
                .Sphere => |_| {
                    switch (shapes[index_b]) {
                        .Sphere => |_| collideSphereSphere(
                            @intCast(index_a), 
                            transforms[index_a], 
                            shapes[index_a],
                            @intCast(index_b), 
                            transforms[index_b], 
                            shapes[index_b],
                            contacts_out
                        ),
                        .Box => |_| collideSphereBox(
                            @intCast(index_a), 
                            transforms[index_a], 
                            shapes[index_a],
                            @intCast(index_b), 
                            transforms[index_b], 
                            shapes[index_b],
                            contacts_out
                        ),
                        else => {},
                    }
                },
                .Box => |_| {
                    switch (shapes[index_b]) {
                        .Sphere => |_| {
                            // sphere-box expects sphere as A, box as B; swap roles
                            const previous_len = contacts_out.items.len;
                            collideSphereBox(
                                @intCast(index_b), 
                                transforms[index_b], 
                                shapes[index_b],
                                @intCast(index_a), 
                                transforms[index_a], 
                                shapes[index_a],
                                contacts_out
                            );

                            // If a contact was added, remap to keep ordering A=index_a (box), B=index_b (sphere)
                            if (contacts_out.items.len > previous_len) {
                                var contact_ref = &contacts_out.items[contacts_out.items.len - 1];

                                // ensure normal points from A(index_a) to B(index_b)
                                contact_ref.normal = contact_ref.normal.negate();

                                // set indices to match (index_a -> index_b)
                                contact_ref.body_a = @intCast(index_a);
                                contact_ref.body_b = @intCast(index_b);

                                // swap world space contact points to preserve A/B semantics
                                const tmp = contact_ref.point_a;
                                contact_ref.point_a = contact_ref.point_b;
                                contact_ref.point_b = tmp;
                            }
                        },
                        .Box => |_| collideBoxBox(
                            @intCast(index_a), 
                            transforms[index_a], 
                            shapes[index_a],
                            @intCast(index_b), 
                            transforms[index_b], 
                            shapes[index_b],
                            manifolds_out
                        ),
                        else => {},
                    }
                },
                else => {},
            }
        }
    }
}

/// Build penetration constraints from contacts
pub fn buildPenetrationConstraints(
    bodies: std.MultiArrayList(body_module.BodyComponents).Slice,
    contacts: []const Contact,
    manifolds: []const ContactManifold,
    constraints_out: *std.ArrayList(contact.PenetrationConstraint)
) void {

    
    // Build constraints from individual contacts
    for (contacts) |contact_entry| {
        const constraint = buildPenetrationConstraint(
            bodies,
            contact_entry.point_a,
            contact_entry.point_b,
            contact_entry.normal,
            contact_entry.body_a,
            contact_entry.body_b,
        );
        constraints_out.appendAssumeCapacity(constraint);
    }
    
    // Build constraints from manifolds
    for (manifolds) |manifold| {
        for (0..manifold.length) |i| {
            const constraint = buildPenetrationConstraint(
                bodies,
                manifold.contact_points_a[i],
                manifold.contact_points_b[i],
                manifold.normal,
                manifold.body_a,
                manifold.body_b,
            );
            constraints_out.appendAssumeCapacity(constraint);
        }
    }
}


/// Helper function to build a single penetration constraint
inline fn buildPenetrationConstraint(
    bodies: std.MultiArrayList(body_module.BodyComponents).Slice,
    contact_point_a: math.Vec3,
    contact_point_b: math.Vec3,
    normal: math.Vec3,
    body_a :u32,
    body_b: u32,
) contact.PenetrationConstraint {
    const motion = bodies.items(.motion);
    const transform = bodies.items(.transform);
    const physics_props = bodies.items(.physics_props);

    const motion_a: MotionComp = motion[body_a];
    const motion_b: MotionComp = motion[body_b];
    const transform_a: TransformComp = transform[body_a];
    const transform_b: TransformComp = transform[body_b];
    const physics_props_a: PhysicsPropsComp = physics_props[body_a];
    const physics_props_b: PhysicsPropsComp = physics_props[body_b];

    const r1 = contact_point_a.sub(&transform_a.position);
    const r2 = contact_point_b.sub(&transform_b.position);

    const inv_mass_a = physics_props_a.inverseMass;
    const inv_mass_b = physics_props_b.inverseMass;

    // Compute tangent directions
    const out_tangent1 = normal.getNormalizePerpendicular();
    const out_tangent2 = normal.cross(&out_tangent1);

    // Compute all cross products
    const r1_cross_n = r1.cross(&normal);
    const r2_cross_n = r2.cross(&normal);
    const r1_cross_t1 = r1.cross(&out_tangent1);
    const r2_cross_t1 = r2.cross(&out_tangent1);
    const r1_cross_t2 = r1.cross(&out_tangent2);
    const r2_cross_t2 = r2.cross(&out_tangent2);

    // Precompute world-space inertia transformations using proper rotation
    const inv_inertia_r1_cross_n = applyDiagonalInertiaWorld(transform_a.orientation, physics_props_a.inverseInertia, r1_cross_n);
    const inv_inertia_r2_cross_n = applyDiagonalInertiaWorld(transform_b.orientation, physics_props_b.inverseInertia, r2_cross_n);
    const inv_inertia_r1_cross_t1 = applyDiagonalInertiaWorld(transform_a.orientation, physics_props_a.inverseInertia, r1_cross_t1);
    const inv_inertia_r2_cross_t1 = applyDiagonalInertiaWorld(transform_b.orientation, physics_props_b.inverseInertia, r2_cross_t1);
    const inv_inertia_r1_cross_t2 = applyDiagonalInertiaWorld(transform_a.orientation, physics_props_a.inverseInertia, r1_cross_t2);
    const inv_inertia_r2_cross_t2 = applyDiagonalInertiaWorld(transform_b.orientation, physics_props_b.inverseInertia, r2_cross_t2);

    // Velocity bias for restitution only
    const relative_vel = motion_b.velocity.sub(&motion_a.velocity);
    const closing_velocity = relative_vel.dot(&normal);

    const restitution = if (closing_velocity < -0.5) @max(physics_props_a.restitution, physics_props_b.restitution) else 0.0;
    const linear_velocity_bias = -restitution * closing_velocity;
    const combined_friction = @sqrt(physics_props_a.friction * physics_props_b.friction);

    // Note: velocity_bias in PenetrationConstraint is Vec3, storing scalar in x component
    return .{
        .r1 = r1,
        .r2 = r2,
        .n = normal,
        .inv_inertia_r1_cross_n = inv_inertia_r1_cross_n,
        .inv_inertia_r2_cross_n = inv_inertia_r2_cross_n,
        .inv_inertia_r1_cross_t1 = inv_inertia_r1_cross_t1,
        .inv_inertia_r2_cross_t1 = inv_inertia_r2_cross_t1,
        .inv_inertia_r1_cross_t2 = inv_inertia_r1_cross_t2,
        .inv_inertia_r2_cross_t2 = inv_inertia_r2_cross_t2,
        .velocity_bias = linear_velocity_bias,
        .accumulated_impulse = 0,
        .accumulated_impulse_tangent1 = 0,
        .accumulated_impulse_tangent2 = 0,
        .inv_mass_a = inv_mass_a,
        .inv_mass_b = inv_mass_b,
        .friction = combined_friction,
        .body_a = body_a,
        .body_b = body_b,
    };
}


/// Apply diagonal inverse inertia in world space to a vector
/// Computes: I_world^(-1) * v = R * I_local^(-1) * R^T * v
inline fn applyDiagonalInertiaWorld(orientation: math.Quat, inv_inertia_local: math.Vec3, v_world: math.Vec3) math.Vec3 {
    // Transform v to local space: R^T * v
    const q_conj = orientation.conjugate();
    const v_local = v_world.mulQuat(&q_conj);

    // Apply diagonal inertia: I_local^(-1) * v_local (component-wise)
    const result_local = inv_inertia_local.mul(&v_local);
    // Transform back to world: R * result_local
    return result_local.mulQuat(&orientation);
}

inline fn effectiveMass(dir: math.Vec3, inv_mass: f32, inv_inertia_world: math.Mat3x3, r_world: math.Vec3) f32 {
    const lever_arm = r_world.cross(&dir);
    const angular_component = inv_inertia_world.mulVec(&lever_arm);
    return inv_mass + lever_arm.dot(&angular_component);
}


// Todo: This position solver needs some work
pub fn solvePosition(
    bodies: std.MultiArrayList(BodyComponents).Slice,
    contacts: []const Contact,
    manifolds: []const ContactManifold,
    iterations: u32
) void {
    _ = iterations;
    const correction_percent: f32 = 0.2;
    const penetration_slop: f32 = 0.02;

    const transform = bodies.items(.transform);
    const physics_props = bodies.items(.physics_props);

    for (contacts) |contact_entry| {
        const inv_mass_a = physics_props[contact_entry.body_a].inverseMass;
        const inv_mass_b = physics_props[contact_entry.body_b].inverseMass;
        const inv_mass_sum = inv_mass_a + inv_mass_b;
        if (inv_mass_sum == 0) continue;

        const correction_magnitude = correction_percent * @max(contact_entry.penetration - penetration_slop, 0.0) / inv_mass_sum;
        const correction = contact_entry.normal.normalize(math.eps_f32).mulScalar(correction_magnitude);

        transform[contact_entry.body_a].position = transform[contact_entry.body_a].position.sub(&correction.mulScalar(inv_mass_a));
        transform[contact_entry.body_b].position = transform[contact_entry.body_b].position.add(&correction.mulScalar(inv_mass_b));
    }

    for (manifolds) |manifold| {
        const inv_mass_a = physics_props[manifold.body_a].inverseMass;
        const inv_mass_b = physics_props[manifold.body_b].inverseMass;
        const inv_mass_sum = inv_mass_a + inv_mass_b;
        if (inv_mass_sum == 0) continue;

        const correction_magnitude = correction_percent * @max(manifold.penetration_depth - penetration_slop, 0.0) / inv_mass_sum;
        const correction = manifold.normal.normalize(math.eps_f32).mulScalar(correction_magnitude);

        transform[manifold.body_a].position = transform[manifold.body_a].position.sub(&correction.mulScalar(inv_mass_a));
        transform[manifold.body_b].position = transform[manifold.body_b].position.add(&correction.mulScalar(inv_mass_b));
    }
}
