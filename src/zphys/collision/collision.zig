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

pub fn solveVelocity(
    bodies: std.MultiArrayList(BodyComponents).Slice,
    contacts: []const Contact,
    manifolds: []const ContactManifold,
    iterations: u32
) void {
    const motion = bodies.items(.motion);
    const transform = bodies.items(.transform);
    const physics_props = bodies.items(.physics_props);
    
    var iteration_index: u32 = 0;
    while (iteration_index < iterations) : (iteration_index += 1) {
        for (contacts) |contact_entry| {
            const contact_normal = contact_entry.normal.normalize(0);
            solveContactPoint(
                &motion[contact_entry.body_a],
                &motion[contact_entry.body_b],
                transform[contact_entry.body_a],
                transform[contact_entry.body_b],
                physics_props[contact_entry.body_a],
                physics_props[contact_entry.body_b],
                contact_normal,
                contact_entry.point_a,
                contact_entry.point_b,
                contact_entry.penetration
            );
        }

        for (manifolds) |manifold| {
            const contact_normal = manifold.normal.normalize(0);
            for (0..manifold.length) |i| {
                solveContactPoint(
                    &motion[manifold.body_a],
                    &motion[manifold.body_b],
                    transform[manifold.body_a],
                    transform[manifold.body_b],
                    physics_props[manifold.body_a],
                    physics_props[manifold.body_b],
                    contact_normal,
                    manifold.contact_points_a[i],
                    manifold.contact_points_b[i],
                    manifold.penetration_depth
                );
            }
        }
    }
}

inline fn solveContactPoint(
    motion_a: *MotionComp,
    motion_b: *MotionComp,
    transform_a: TransformComp,
    transform_b: TransformComp,
    physics_props_a: PhysicsPropsComp,
    physics_props_b: PhysicsPropsComp,
    contact_normal: math.Vec3,
    point_world_a: math.Vec3,
    point_world_b: math.Vec3,
    penetration: f32
) void {
    const penetration_slop: f32 = 0.003;

    // Compute lever arms from body centers to contact points (world space)
    const r_a_world = point_world_a.sub(&transform_a.position);
    const r_b_world = point_world_b.sub(&transform_b.position);
    const vel_a = motion_a.velocity.add(&motion_a.angularVelocity.cross(&r_a_world));
    const vel_b = motion_b.velocity.add(&motion_b.angularVelocity.cross(&r_b_world));
    const relative_velocity = vel_b.sub(&vel_a);
    const velocity_along_normal = relative_velocity.dot(&contact_normal);

    const corrected_penetration = @max(penetration - penetration_slop, 0.0);
    // skip when objects are moving apart and barely touching
    if (velocity_along_normal > 0 and corrected_penetration <= 0) return;

    // Restitution only on closing velocity
    const restitution = if (velocity_along_normal < -0.5) @max(physics_props_a.restitution, physics_props_b.restitution) else 0.0;

    const inv_mass_a = physics_props_a.inverseMass;
    const inv_mass_b = physics_props_b.inverseMass;
    const inv_inertia_world_a = computeInverseInertiaWorld(transform_a, physics_props_a);
    const inv_inertia_world_b = computeInverseInertiaWorld(transform_b, physics_props_b);

    // Normal impulse (include restitution + bias)
    const k_a_n = effectiveMass(contact_normal, inv_mass_a, inv_inertia_world_a, r_a_world);
    const k_b_n = effectiveMass(contact_normal, inv_mass_b, inv_inertia_world_b, r_b_world);
    const k_n = k_a_n + k_b_n;
    if (k_n <= 0) return;
    var normal_impulse_magnitude = (-(1.0 + restitution) * velocity_along_normal) / k_n;
    if (normal_impulse_magnitude < 0) normal_impulse_magnitude = 0;

    if (normal_impulse_magnitude > 0) {
        const normal_impulse = contact_normal.mulScalar(normal_impulse_magnitude);
        motion_a.velocity = motion_a.velocity.sub(&normal_impulse.mulScalar(inv_mass_a));
        motion_b.velocity = motion_b.velocity.add(&normal_impulse.mulScalar(inv_mass_b));

        const angular_impulse_a = r_a_world.cross(&normal_impulse);
        const delta_omega_a = inv_inertia_world_a.mulVec(&angular_impulse_a);
        motion_a.angularVelocity = motion_a.angularVelocity.sub(&delta_omega_a);

        const angular_impulse_b = r_b_world.cross(&normal_impulse);
        const delta_omega_b = inv_inertia_world_b.mulVec(&angular_impulse_b);
        motion_b.angularVelocity = motion_b.angularVelocity.add(&delta_omega_b);
    }

    // Friction (Coulomb, clamped by mu * |jn|)
    var vel_tan = relative_velocity.sub(&contact_normal.mulScalar(relative_velocity.dot(&contact_normal)));
    const tangent_len2 = vel_tan.len2();
    if (tangent_len2 <= 1e-12) return;
    vel_tan = vel_tan.normalize(math.eps_f32);

    // Effective mass in tangent direction
    const k_a_t = effectiveMass(vel_tan, inv_mass_a, inv_inertia_world_a, r_a_world);
    const k_b_t = effectiveMass(vel_tan, inv_mass_b, inv_inertia_world_b, r_b_world);
    const k_t = k_a_t + k_b_t;
    const tangential_impulse_magnitude = -(relative_velocity.dot(&vel_tan)) / k_t;
    const friction_coefficient = std.math.sqrt(@max(physics_props_a.friction, 0) * @max(physics_props_b.friction, 0));

    const max_friction = friction_coefficient * normal_impulse_magnitude;
    var clamped_tangential_impulse = tangential_impulse_magnitude;
    if (clamped_tangential_impulse > max_friction) clamped_tangential_impulse = max_friction;
    if (clamped_tangential_impulse < -max_friction) clamped_tangential_impulse = -max_friction;

    const tangential_impulse = vel_tan.mulScalar(clamped_tangential_impulse);

    // Linear friction impulses
    motion_a.velocity = motion_a.velocity.sub(&tangential_impulse.mulScalar(inv_mass_a));
    motion_b.velocity = motion_b.velocity.add(&tangential_impulse.mulScalar(inv_mass_b));

    // Angular friction impulses
    const angular_impulse_a_t = r_a_world.cross(&tangential_impulse);
    const delta_omega_a_t = inv_inertia_world_a.mulVec(&angular_impulse_a_t);
    motion_a.angularVelocity = motion_a.angularVelocity.sub(&delta_omega_a_t);
    const angular_impulse_b_t = r_b_world.cross(&tangential_impulse);
    const delta_omega_b_t = inv_inertia_world_b.mulVec(&angular_impulse_b_t);
    motion_b.angularVelocity = motion_b.angularVelocity.add(&delta_omega_b_t);
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
    const n = normal.normalize(0);

    const inv_mass_a = physics_props_a.inverseMass;
    const inv_mass_b = physics_props_b.inverseMass;
    const inv_inertia_a = computeInverseInertiaWorld(transform_a, physics_props_a);
    const inv_inertia_b = computeInverseInertiaWorld(transform_b, physics_props_b);

    // Precompute inverse inertia × (r × n)
    const r1_cross_n = r1.cross(&n);
    const r2_cross_n = r2.cross(&n);
    const inv_inertia_r1_cross_n = inv_inertia_a.mulVec(&r1_cross_n);
    const inv_inertia_r2_cross_n = inv_inertia_b.mulVec(&r2_cross_n);

    // Compute effective mass
    const k1 = inv_mass_a + inv_mass_b;
    const k2 = inv_inertia_r1_cross_n.dot(&r1_cross_n);
    const k3 = inv_inertia_r2_cross_n.dot(&r2_cross_n);
    const inverse_effective_mass_val = k1 + k2 + k3;

    // Velocity bias for restitution only
    const vel_a = motion_a.velocity.add(&motion_a.angularVelocity.cross(&r1));
    const vel_b = motion_b.velocity.add(&motion_b.angularVelocity.cross(&r2));
    const relative_vel = vel_b.sub(&vel_a);
    const closing_velocity = relative_vel.dot(&n);

    const restitution = if (closing_velocity < -0.5) @max(physics_props_a.restitution, physics_props_b.restitution) else 0.0;

    const velocity_bias_val = -restitution * closing_velocity;

    // Note: velocity_bias in PenetrationConstraint is Vec3, storing scalar in x component
    return .{
        .r1 = r1,
        .r2 = r2,
        .n = n,
        .velocity_bias = velocity_bias_val,
        .invert_inertia_n_x_r1 = inv_inertia_r1_cross_n,
        .invert_inertia_n_x_r2 = inv_inertia_r2_cross_n,
        .inverse_effective_mass = if (inverse_effective_mass_val > 0) 1.0 / inverse_effective_mass_val else 0,
        .accumulated_impulse = 0,
        .inv_mass_a = inv_mass_a,
        .inv_mass_b = inv_mass_b,
        .body_a = body_a,
        .body_b = body_b,
    };
}

inline fn computeInverseInertiaWorld(transform: TransformComp, physics_props: PhysicsPropsComp) math.Mat3x3 {
    // Build world-space inverse inertia once from local and orientation
    const q = transform.orientation.normalize();
    const rot4 = math.Mat4x4.rotateByQuaternion(q);
    const r0 = rot4.row(0);
    const r1 = rot4.row(1);
    const r2 = rot4.row(2);
    const rot3 = math.Mat3x3.init(
        &math.vec3(r0.x(), r0.y(), r0.z()),
        &math.vec3(r1.x(), r1.y(), r1.z()),
        &math.vec3(r2.x(), r2.y(), r2.z()),
    );
    const rot3_t = rot3.transpose();
    return rot3.mul(&physics_props.inverseInertia).mul(&rot3_t);
}

inline fn effectiveMass(dir: math.Vec3, inv_mass: f32, inv_inertia_world: math.Mat3x3, r_world: math.Vec3) f32 {
    const lever_arm = r_world.cross(&dir);
    const angular_component = inv_inertia_world.mulVec(&lever_arm);
    return inv_mass + lever_arm.dot(&angular_component);
}

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
        const correction = contact_entry.normal.normalize(0).mulScalar(correction_magnitude);

        transform[contact_entry.body_a].position = transform[contact_entry.body_a].position.sub(&correction.mulScalar(inv_mass_a));
        transform[contact_entry.body_b].position = transform[contact_entry.body_b].position.add(&correction.mulScalar(inv_mass_b));
    }

    for (manifolds) |manifold| {
        const inv_mass_a = physics_props[manifold.body_a].inverseMass;
        const inv_mass_b = physics_props[manifold.body_b].inverseMass;
        const inv_mass_sum = inv_mass_a + inv_mass_b;
        if (inv_mass_sum == 0) continue;

        const correction_magnitude = correction_percent * @max(manifold.penetration_depth - penetration_slop, 0.0) / inv_mass_sum;
        const correction = manifold.normal.normalize(0).mulScalar(correction_magnitude);

        transform[manifold.body_a].position = transform[manifold.body_a].position.sub(&correction.mulScalar(inv_mass_a));
        transform[manifold.body_b].position = transform[manifold.body_b].position.add(&correction.mulScalar(inv_mass_b));
    }
}
