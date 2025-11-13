const std = @import("std");
const math = @import("math");
const TransformComp = @import("../body.zig").TransformComp;
const Shape = @import("shape.zig").Shape;
const contact = @import("contact.zig");

// Expects: A is Sphere, B is Box
pub fn collideSphereBox(
    a_id: u32, 
    transform_a: TransformComp, 
    shape_a: Shape,
    b_id: u32, 
    transform_b: TransformComp, 
    shape_b: Shape,
    out: *std.ArrayList(contact.Contact)
) void {
    const sphere = shape_a.Sphere;
    const box = shape_b.Box;

    const closest = closestPointOnOBB(transform_a.position, transform_b.position, transform_b.orientation, box.half_extents);
    const vector_box_to_sphere = closest.sub(&transform_a.position);
    const distance_squared = vector_box_to_sphere.len2();

    if (distance_squared > sphere.radius * sphere.radius) return;

    var normal: math.Vec3 = undefined;
    const distance = std.math.sqrt(distance_squared);
    if (distance > 1e-6) {
        normal = vector_box_to_sphere.mulScalar(1.0 / distance); // from box->sphere
    } else {
        // Choose a reasonable normal (up)
        normal = math.vec3(0, 1, 0);
    }

    const penetration = sphere.radius - distance;

    // Contact points stored in world space
    const point_a = transform_a.position.add(&normal.negate().mulScalar(sphere.radius));
    const point_b = closest;

    out.appendAssumeCapacity(.{
        .body_a = a_id,
        .body_b = b_id,
        .normal = normal, // from box to sphere
        .point_a = point_a,
        .point_b = point_b,
        .penetration = penetration,
    });
}

fn closestPointOnOBB(point: math.Vec3, center: math.Vec3, orientation: math.Quat, half_extents: math.Vec3) math.Vec3 {
    const p_local = point.sub(&center).mulQuat(&orientation.conjugate());
    const clamped = math.vec3(
        std.math.clamp(p_local.x(), -half_extents.x(), half_extents.x()),
        std.math.clamp(p_local.y(), -half_extents.y(), half_extents.y()),
        std.math.clamp(p_local.z(), -half_extents.z(), half_extents.z()),
    );
    // Back to world
    return center.add(&clamped.mulQuat(&orientation));
}
