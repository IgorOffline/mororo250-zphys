const std = @import("std");
const math = @import("math");
const Body = @import("../body.zig").Body;
const contact = @import("contact.zig");

pub fn collideSphereSphere(a_id: u32, sphereBodyA: *const Body, b_id: u32, sphereBodyB: *const Body, out: *std.ArrayList(contact.Contact)) void {
    const sphere_a = sphereBodyA.shape.Sphere;
    const sphere_b = sphereBodyB.shape.Sphere;

    const vector_a_to_b = sphereBodyB.position.sub(&sphereBodyA.position);
    const distance_squared = vector_a_to_b.len2();
    const combined_radius = sphere_a.radius + sphere_b.radius;

    if (distance_squared > combined_radius * combined_radius) return; // no contact

    const distance = std.math.sqrt(distance_squared);
    var normal: math.Vec3 = undefined;
    if (distance > 1e-6) {
        normal = vector_a_to_b.mulScalar(1.0 / distance);
    } else {
        normal = math.vec3(0, 1, 0);
    }

    const penetration = combined_radius - distance;
    // Approx contact point as point on surface of A along normal
    // Todo: Calculate two differnt points for each of the objects?
    const point = sphereBodyA.position.add(&normal.mulScalar(sphere_a.radius - penetration * 0.5));

    // I am rotating here to rotate back after
    // Todo: Save as world space and calculate anything need in local space later
    // This way we avoid half of the calcualtions at least
    const inv_q_a = sphereBodyA.orientation.conjugate();
    const inv_q_b = sphereBodyB.orientation.conjugate();
    const point_local_a = point.sub(&sphereBodyA.position).mulQuat(&inv_q_a);
    const point_local_b = point.sub(&sphereBodyB.position).mulQuat(&inv_q_b);

    out.appendAssumeCapacity( .{
        .body_a = a_id,
        .body_b = b_id,
        .normal = normal,
        .point_local_a = point_local_a,
        .point_local_b = point_local_b,
        .penetration = penetration,
    });
}
