const std = @import("std");
const math = @import("math");
const TransformComp = @import("../body.zig").TransformComp;
const Shape = @import("shape.zig").Shape;
const contact = @import("contact.zig");
const gjk = @import("gjk.zig");
const epa = @import("epa.zig");
const sat = @import("sat.zig");
const manifold_between_two_faces = @import("manifold_between_two_faces.zig");

pub fn collideBoxBox(
    a_id: u32, 
    transform_a: TransformComp, 
    shape_a: Shape,
    b_id: u32, 
    transform_b: TransformComp, 
    shape_b: Shape,
    out: *std.ArrayList(contact.ContactManifold)
) void {
    const box_a = shape_a.Box;
    const box_b = shape_b.Box;

    // GJK for detection
    const gjk_shape_a = gjk.GjkBox{ .center = transform_a.position, .orientation = transform_a.orientation, .half_extents = box_a.half_extents };
    const gjk_shape_b = gjk.GjkBox{ .center = transform_b.position, .orientation = transform_b.orientation, .half_extents = box_b.half_extents };

    // CSO (A − B) and support-point buffers:
    // - For box–box, the CSO can have up to 16 vertices in general position.
    // - EPA grows the Minkowski simplex by appending support points; we must also
    //   store the matching A/B support points for each Minkowski vertex.
    // - Therefore, all three arrays below must have the same capacity (16).
    var minkowski_points: [16]math.Vec3 = undefined;
    var shape_a_points: [16]math.Vec3 = undefined;
    var shape_b_points: [16]math.Vec3 = undefined;
    // Note: With V = 16, EPA's worst-case face count is F ≤ 2V − 4 = 28 (see epa.zig).
    const simplex_arrays: [3][]math.Vec3 = .{ minkowski_points[0..], shape_a_points[0..], shape_b_points[0..] };

    const intersects = gjk.gjkIntersect(simplex_arrays, gjk_shape_a, gjk_shape_b);
    if (!intersects) return;

    const epa_result = epa.epa(simplex_arrays, gjk_shape_a, gjk_shape_b);

    var penetration_axis = epa_result.penetration_axis;
    const delta_centers = transform_b.position.sub(&transform_a.position);
    if (delta_centers.dot(&penetration_axis) < 0) {
        penetration_axis = penetration_axis.negate();
    }

    const face_a = gjk_shape_a.getSupportFace(penetration_axis);
    const face_b = gjk_shape_b.getSupportFace(penetration_axis.negate());

    const max_length = face_a.len + face_b.len;
    var face_a_contact_points: [max_length]math.Vec3 = undefined;
    var face_b_contact_points: [max_length]math.Vec3 = undefined;

    const manifold_size = manifold_between_two_faces.manifoldBetweenTwoFaces(
        face_a.len + face_b.len,
        &face_a,
        &face_b,
        penetration_axis,
        &face_a_contact_points,
        &face_b_contact_points,
    ) catch {
        // Fallback: store world space points directly
        out.appendAssumeCapacity(.{
            .body_a = a_id,
            .body_b = b_id,
            .normal = penetration_axis.normalize(math.eps_f32),
            .penetration_depth = epa_result.penetration_depth,
            .length = 1,
            .contact_points_a = .{ epa_result.collision_point_a, undefined, undefined, undefined },
            .contact_points_b = .{ epa_result.collision_point_b, undefined, undefined, undefined },
        });
        return;
    };

    out.appendAssumeCapacity(.{
        .body_a = a_id,
        .body_b = b_id,
        .normal = penetration_axis.normalize(math.eps_f32),
        .penetration_depth = epa_result.penetration_depth,
        .contact_points_a = undefined,
        .contact_points_b = undefined,
        .length = undefined,
    });

    var manifold: *contact.ContactManifold = &out.items[out.items.len - 1];
    var contact_points_a: []math.Vec3 = &manifold.contact_points_a;
    var contact_points_b: []math.Vec3 = &manifold.contact_points_b;

    if (manifold_size > 4) {
        manifold_between_two_faces.pruneContactPoints(
            max_length,
            penetration_axis,
            face_a_contact_points[0..manifold_size],
            face_b_contact_points[0..manifold_size],
            &contact_points_a,
            &contact_points_b,
        );
        manifold.length = @intCast(contact_points_a.len);
    } else {
        // Copy directly when <= 4 points
        for (0..manifold_size) |i| {
            manifold.contact_points_a[i] = face_a_contact_points[i];
            manifold.contact_points_b[i] = face_b_contact_points[i];
        }
        manifold.length = @intCast(manifold_size);
    }
}
