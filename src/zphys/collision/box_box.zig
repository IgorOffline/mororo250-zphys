const std = @import("std");
const math = @import("math");
const Body = @import("../body.zig").Body;
const contact = @import("contact.zig");
const gjk = @import("gjk.zig");
const epa = @import("epa.zig");
const sat = @import("sat.zig");

pub fn collideBoxBox(a_id: u32, body_a: *const Body, b_id: u32, body_b: *const Body, out: *std.ArrayList(contact.Contact)) void {
    const box_a = body_a.shape.Box;
    const box_b = body_b.shape.Box;

    // GJK for detection
    const shape_a = gjk.GjkBox{ .center = body_a.position, .orientation = body_a.orientation, .half_extents = box_a.half_extents };
    const shape_b = gjk.GjkBox{ .center = body_b.position, .orientation = body_b.orientation, .half_extents = box_b.half_extents };
    var minkowski_points: [16]math.Vec3 = undefined;
    var shape_a_points: [8]math.Vec3 = undefined;
    var shape_b_points: [8]math.Vec3 = undefined;
    const simplex_arrays: [3][]math.Vec3 = .{ minkowski_points[0..], shape_a_points[0..], shape_b_points[0..] };

    const intersects = gjk.gjkIntersect(simplex_arrays, shape_a, shape_b);
    if (!intersects) return;

    const epa_result = epa.epa(simplex_arrays, shape_a, shape_b);

    const inv_q_a = body_a.orientation.inverse();
    const inv_q_b = body_b.orientation.inverse();
    const point_local_a = epa_result.collision_point_a.sub(&body_a.position).mulQuat(&inv_q_a);
    const point_local_b = epa_result.collision_point_b.sub(&body_b.position).mulQuat(&inv_q_b);

    var n = epa_result.normal;
    const delta_centers = body_b.position.sub(&body_a.position);
    if (delta_centers.dot(&n) < 0) {
        n = n.negate();
    }

    out.appendAssumeCapacity( .{
        .body_a = a_id,
        .body_b = b_id,
        .normal = n,
        .point_local_a = point_local_a,
        .point_local_b = point_local_b,
        .penetration = epa_result.penetration_depth,
    });
}
