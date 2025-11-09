const std = @import("std");
const math = @import("math");
const Body = @import("../body.zig").Body;

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

pub const gjkBoxesIntersect = gjk.gjkBoxesIntersect;
pub const satBoxBoxContact = sat.satBoxBoxContact;
pub const ContactManifold = contact.ContactManifold;

// Todo: Add BroadPhase collision check in here
pub fn generateContacts(bodies: []const Body, contacts_out: *std.ArrayList(Contact), manifolds_out: *std.ArrayList(contact.ContactManifold)) void {
    var index_a: usize = 0;
    while (index_a < bodies.len) : (index_a += 1) {
        var index_b: usize = index_a + 1;
        while (index_b < bodies.len) : (index_b += 1) {
            const body_a = &bodies[index_a];
            const body_b = &bodies[index_b];

            if (body_a.inverseMass == 0 and body_b.inverseMass == 0) continue;

            switch (body_a.shape) {
                .Sphere => |_| {
                    switch (body_b.shape) {
                        .Sphere => |_| collideSphereSphere(@intCast(index_a), body_a, @intCast(index_b), body_b, contacts_out),
                        .Box => |_| collideSphereBox(@intCast(index_a), body_a, @intCast(index_b), body_b, contacts_out),
                        else => {},
                    }
                },
                .Box => |_| {
                    switch (body_b.shape) {
                        .Sphere => |_| {
                            // sphere-box expects sphere as A, box as B; swap roles
                            const previous_len = contacts_out.items.len;
                            collideSphereBox(@intCast(index_b), body_b, @intCast(index_a), body_a, contacts_out);

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
                        .Box => |_| collideBoxBox(@intCast(index_a), body_a, @intCast(index_b), body_b, manifolds_out),
                        else => {},
                    }
                },
                else => {},
            }
        }
    }
}

pub fn solveVelocity(bodies: []Body, contacts: []const Contact, manifolds: []const ContactManifold, iterations: u32) void {
    var iteration_index: u32 = 0;
    while (iteration_index < iterations) : (iteration_index += 1) {
        for (contacts) |contact_entry| {
            const body_a = &bodies[contact_entry.body_a];
            const body_b = &bodies[contact_entry.body_b];
            const contact_normal = contact_entry.normal.normalize(0);

            solveContactPoint(body_a, body_b, contact_normal, contact_entry.point_a, contact_entry.point_b, contact_entry.penetration);
        }

        for (manifolds) |manifold| {
            const body_a = &bodies[manifold.body_a];
            const body_b = &bodies[manifold.body_b];
            const contact_normal = manifold.normal.normalize(0);

            for (0..manifold.length) |i| {
                solveContactPoint(body_a, body_b, contact_normal, manifold.contact_points_a[i], manifold.contact_points_b[i], manifold.penetration_depth);
            }
        }
    }
}

inline fn solveContactPoint(body_a: *Body, body_b: *Body, contact_normal: math.Vec3, point_world_a: math.Vec3, point_world_b: math.Vec3, penetration: f32) void {
    const penetration_slop: f32 = 0.003;

    // Compute lever arms from body centers to contact points (world space)
    const r_a_world = point_world_a.sub(&body_a.position);
    const r_b_world = point_world_b.sub(&body_b.position);
    const vel_a = body_a.velocity.add(&body_a.angularVelocity.cross(&r_a_world));
    const vel_b = body_b.velocity.add(&body_b.angularVelocity.cross(&r_b_world));
    const relative_velocity = vel_b.sub(&vel_a);
    const velocity_along_normal = relative_velocity.dot(&contact_normal);

    const corrected_penetration = @max(penetration - penetration_slop, 0.0);
    // skip when objects are moving a part and barely touching
    if (velocity_along_normal > 0 and corrected_penetration <= 0) return;

    // Restitution only on closing velocity
    const restitution = if (velocity_along_normal < -0.5) @max(body_a.restitution, body_b.restitution) else 0.0;

    const inv_mass_a = body_a.inverseMass;
    const inv_mass_b = body_b.inverseMass;
    const inv_inertia_world_a = computeInverseInertiaWorld(body_a);
    const inv_inertia_world_b = computeInverseInertiaWorld(body_b);

    // Normal impulse (include restitution + bias)
    const k_a_n = effectiveMass(contact_normal, inv_mass_a, inv_inertia_world_a, r_a_world);
    const k_b_n = effectiveMass(contact_normal, inv_mass_b, inv_inertia_world_b, r_b_world);
    const k_n = k_a_n + k_b_n;
    if (k_n <= 0) return;
    var normal_impulse_magnitude = (-(1.0 + restitution) * velocity_along_normal) / k_n;
    if (normal_impulse_magnitude < 0) normal_impulse_magnitude = 0;

    if (normal_impulse_magnitude > 0) {
        const normal_impulse = contact_normal.mulScalar(normal_impulse_magnitude);
        body_a.velocity = body_a.velocity.sub(&normal_impulse.mulScalar(inv_mass_a));
        body_b.velocity = body_b.velocity.add(&normal_impulse.mulScalar(inv_mass_b));

        const angular_impulse_a = r_a_world.cross(&normal_impulse);
        const delta_omega_a = inv_inertia_world_a.mulVec(&angular_impulse_a);
        body_a.angularVelocity = body_a.angularVelocity.sub(&delta_omega_a);

        const angular_impulse_b = r_b_world.cross(&normal_impulse);
        const delta_omega_b = inv_inertia_world_b.mulVec(&angular_impulse_b);
        body_b.angularVelocity = body_b.angularVelocity.add(&delta_omega_b);
    }

    // Friction (Coulomb, clamped by mu * |jn|)
    var vel_tan = relative_velocity.sub(&contact_normal.mulScalar(relative_velocity.dot(&contact_normal)));
    const tangent_len2 = vel_tan.len2();
    if (tangent_len2 <= 1e-12) return;
    vel_tan = vel_tan.normalize(math.eps_f32);

    // Effective mass in tangent direction
    const k_a_t = effectiveMass(vel_tan, inv_mass_a, body_a.inverseInertia, r_a_world);
    const k_b_t = effectiveMass(vel_tan, inv_mass_b, body_b.inverseInertia, r_b_world);
    const k_t = k_a_t + k_b_t;
    const tangential_impulse_magnitude = -(relative_velocity.dot(&vel_tan)) / k_t;
    const friction_coefficient = std.math.sqrt(@max(body_a.friction, 0) * @max(body_b.friction, 0));

    const max_friction = friction_coefficient * normal_impulse_magnitude;
    var clamped_tangential_impulse = tangential_impulse_magnitude;
    if (clamped_tangential_impulse > max_friction) clamped_tangential_impulse = max_friction;
    if (clamped_tangential_impulse < -max_friction) clamped_tangential_impulse = -max_friction;

    const tangential_impulse = vel_tan.mulScalar(clamped_tangential_impulse);

    // Linear friction impulses
    body_a.velocity = body_a.velocity.sub(&tangential_impulse.mulScalar(inv_mass_a));
    body_b.velocity = body_b.velocity.add(&tangential_impulse.mulScalar(inv_mass_b));

    // Angular friction impulses
    const angular_impulse_a_t = r_a_world.cross(&tangential_impulse);
    const delta_omega_a_t = body_a.inverseInertia.mulVec(&angular_impulse_a_t);
    body_a.angularVelocity = body_a.angularVelocity.sub(&delta_omega_a_t);
    const angular_impulse_b_t = r_b_world.cross(&tangential_impulse);
    const delta_omega_b_t = body_b.inverseInertia.mulVec(&angular_impulse_b_t);
    body_b.angularVelocity = body_b.angularVelocity.add(&delta_omega_b_t);
}

inline fn computeInverseInertiaWorld(body: *const Body) math.Mat3x3 {
    // Build world-space inverse inertia once from local and orientation
    const q = body.orientation.normalize();
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
    return rot3.mul(&body.inverseInertia).mul(&rot3_t);
}

inline fn effectiveMass(dir: math.Vec3, inv_mass: f32, inv_inertia_world: math.Mat3x3, r_world: math.Vec3) f32 {
    const lever_arm = r_world.cross(&dir);
    const angular_component = inv_inertia_world.mulVec(&lever_arm);
    return inv_mass + lever_arm.dot(&angular_component);
}

pub fn solvePosition(bodies: []Body, contacts: []const Contact, manifolds: []const ContactManifold) void {
    const correction_percent: f32 = 0.2;
    const penetration_slop: f32 = 0.02;

    for (contacts) |contact_entry| {
        const body_a = &bodies[contact_entry.body_a];
        const body_b = &bodies[contact_entry.body_b];

        const inv_mass_a = body_a.inverseMass;
        const inv_mass_b = body_b.inverseMass;
        const inv_mass_sum = inv_mass_a + inv_mass_b;
        if (inv_mass_sum == 0) continue;

        const correction_magnitude = correction_percent * @max(contact_entry.penetration - penetration_slop, 0.0) / inv_mass_sum;
        const correction = contact_entry.normal.normalize(0).mulScalar(correction_magnitude);

        body_a.position = body_a.position.sub(&correction.mulScalar(inv_mass_a));
        body_b.position = body_b.position.add(&correction.mulScalar(inv_mass_b));
    }

    for (manifolds) |manifold| {
        const body_a = &bodies[manifold.body_a];
        const body_b = &bodies[manifold.body_b];

        const inv_mass_a = body_a.inverseMass;
        const inv_mass_b = body_b.inverseMass;
        const inv_mass_sum = inv_mass_a + inv_mass_b;
        if (inv_mass_sum == 0) continue;

        const correction_magnitude = correction_percent * @max(manifold.penetration_depth - penetration_slop, 0.0) / inv_mass_sum;
        const correction = manifold.normal.normalize(0).mulScalar(correction_magnitude);

        body_a.position = body_a.position.sub(&correction.mulScalar(inv_mass_a));
        body_b.position = body_b.position.add(&correction.mulScalar(inv_mass_b));
    }
}
