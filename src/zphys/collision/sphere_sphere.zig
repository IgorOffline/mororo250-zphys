const std = @import("std");
const math = @import("math");
const TransformComp = @import("../body.zig").TransformComp;
const Shape = @import("shape.zig").Shape;
const contact = @import("contact.zig");

pub fn collideSphereSphere(
    a_id: u32, 
    transform_a: TransformComp, 
    shape_a: Shape,
    b_id: u32, 
    transform_b: TransformComp, 
    shape_b: Shape,
    out: *std.ArrayList(contact.Contact)
) void {
    const sphere_a = shape_a.Sphere;
    const sphere_b = shape_b.Sphere;

    const vector_a_to_b = transform_b.position.sub(&transform_a.position);
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

    // Contact points stored in world space
    const point_a = transform_a.position.add(&normal.mulScalar(sphere_a.radius));
    const point_b = transform_b.position.sub(&normal.mulScalar(sphere_b.radius));

    out.appendAssumeCapacity(.{
        .body_a = a_id,
        .body_b = b_id,
        .normal = normal,
        .point_a = point_a,
        .point_b = point_b,
        .penetration = penetration,
    });
}
