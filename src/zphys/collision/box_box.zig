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
    var minkowski_points: [8]math.Vec3 = undefined;
    var shape_a_points: [8]math.Vec3 = undefined;
    var shape_b_points: [8]math.Vec3 = undefined;
    const simplex_arrays: [3][]math.Vec3 = .{ minkowski_points[0..], shape_a_points[0..], shape_b_points[0..] };

    const intersects = gjk.gjkIntersect(simplex_arrays, shape_a, shape_b);
    if (!intersects) return;

    const epa_result = epa.epa(simplex_arrays, shape_a, shape_b);

    const friction = std.math.sqrt(@max(body_a.friction, 0) * @max(body_b.friction, 0));
    const restitution = @max(body_a.restitution, body_b.restitution);

    out.appendAssumeCapacity( .{
        .body_a = a_id,
        .body_b = b_id,
        .normal = epa_result.normal,
        .contact_point_a = epa_result.collision_point_a,
        .contact_point_b = epa_result.collision_point_b,
        .penetration = epa_result.penetration_depth,
        .friction = friction,
        .restitution = restitution,
    });
}
