const std = @import("std");
const math = @import("math");
const clipPoly = @import("sutherland_hodgman.zig");

// We use sutherland hodgman algorithm to clip the faces and find the manifold similar to jolt solution
// Physx does a projection to 2d before prune -> Todo: Compare both solutions later.
pub fn manifoldBetweenTwoFaces(comptime max_length: usize, face_a: []const math.Vec3, face_b: []const math.Vec3, penetration_axis: math.Vec3, out_face_a_contact_points: *[max_length]math.Vec3, out_face_b_contact_points: *[max_length]math.Vec3) !usize {
    const plane_origin = face_a[0];
    const first_edge = face_a[1].sub(&plane_origin);
    const second_edge = face_a[2].sub(&plane_origin);
    const plane_normal = first_edge.cross(&second_edge).normalize(math.eps(f32));
    const penetration_axis_dot_plane_normal = penetration_axis.dot(&plane_normal);

    // If penetration axis and plane normal are perpendicular, fall back to the contact points
    // penetration_axis_dot_pane_normal will be use as the denominator when we calculate the distance from the clipped
    // face to the the project point in the clipping face. This check is necessary for division safety
    if (penetration_axis_dot_plane_normal == 0.0) return error.penetration_perp_plane_normal;

    const clipped_face = clipPoly.clipPolyPoly(max_length, face_b, face_a, penetration_axis.normalize(math.eps_f32), out_face_b_contact_points);
    const penetration_axis_len = std.math.sqrt(penetration_axis.len2());

    // projection step
    // After clipping, the new vertices are on Face 2's surface, but we need contact points on both Face 1 and Face 2
    // to resolve the collision. To solve this problem The projection step finds where each clipped vertex would land
    // on face 1's plane if moved along the penetration direction
    var i: usize = 0;
    var manifold_count: usize = 0;
    while (i < clipped_face.len) : (i += 1) {
        const vertex2 = clipped_face[i];
        // Project clipped face back onto the plane of face 1, we do this by solving:
        // p1 = p2 + distance * penetration_axis / |penetration_axis|
        // (p0 - plane_origin) . plane_normal = 0
        // This gives us:
        // distance = -|penetration_axis| * (p2 - plane_origin) . plane_normal / penetration_axis . plane_normal
        const distance = vertex2.sub(&plane_origin).dot(&plane_normal) / penetration_axis_dot_plane_normal; // note left out -|penetration_axis| term
        const manifold_tolerance = 0.02;
        if (distance * penetration_axis_len < manifold_tolerance) {
            out_face_b_contact_points[manifold_count] = out_face_b_contact_points[i];
            out_face_a_contact_points[manifold_count] = vertex2.sub(&penetration_axis.mulScalar(distance));
            manifold_count += 1;
        }
    }

    // If not contact point was found fall back to old contact point
    if (manifold_count == 0) {
        return error.nocontactpointfound;
    }
    return manifold_count;
}

// This prune solution is the same solution used by jolt physics engine.
// todo: Is there a better heuristics? we should research latter for better possibly cheaper implementations
pub fn pruneContactPoints(comptime max_length: usize, penetration_axis: math.Vec3, contact_points1: []math.Vec3, contact_points2: []math.Vec3, out_contact_point1: *[]math.Vec3, out_contact_point2: *[]math.Vec3) void {
    const min_dist_sq = 1.0e-6;

    var penetration_depth_sq: [max_length]f32 = undefined;
    var projected: [max_length]math.Vec3 = undefined;
    var i: usize = 0;
    while (i < contact_points1.len) : (i += 1) {
        const vertex1 = contact_points1[i];
        projected[i] = vertex1.sub(&penetration_axis.mulScalar(vertex1.dot(&penetration_axis)));
        const vertex2 = contact_points2[i];
        penetration_depth_sq[i] = @max(min_dist_sq, (vertex2.sub(&vertex1).len2()));
    }

    // Use heuristic to find point that is furthest away and has the deepest penetration depth
    var point1_index: usize = 0;
    var heuristic_val = @max(min_dist_sq, projected[0].len2() * penetration_depth_sq[0]);
    i = 1;
    while (i < contact_points1.len) : (i = i + 1) {
        const new_heuristics = @max(min_dist_sq, projected[i].len2() * penetration_depth_sq[i]);
        if (new_heuristics > heuristic_val) {
            heuristic_val = new_heuristics;
            point1_index = i;
        }
    }
    const point1 = projected[point1_index];

    // Combine be far from the first point in the heuristics to look for the second point
    var point2_index: usize = std.math.maxInt(usize);
    heuristic_val = -std.math.floatMax(f32);
    i = 0;
    while (i < contact_points1.len) : (i += 1) {
        if (i == point1_index)
            continue;

        const new_heuristics = @max(min_dist_sq, (projected[i].len2() - point1.len2()) * penetration_depth_sq[i]);
        if (new_heuristics > heuristic_val) {
            heuristic_val = new_heuristics;
            point2_index = i;
        }
    }
    const point2 = projected[point2_index];

    // find furthest points on both sides of the line segment in order to maximize the area
    var point3_index: usize = std.math.maxInt(usize);
    var point4_index: usize = std.math.maxInt(usize);
    var min_dis: f32 = 0;
    var max_dis: f32 = 0;
    i = 0;
    const perp = point2.sub(&point1).cross(&penetration_axis);
    while (i < contact_points1.len) : (i += 1) {
        if (i != point1_index and i != point2_index) {
            const dis = perp.dot(&projected[i].sub(&point1));
            if (dis < min_dis) {
                min_dis = dis;
                point3_index = i;
                continue;
            }
            if (dis > max_dis) {
                max_dis = dis;
                point4_index = i;
            }
        }
    }

    var len: usize = 2;
    out_contact_point1.*[0] = contact_points1[point1_index];
    out_contact_point2.*[0] = contact_points2[point1_index];
    if (point3_index != std.math.maxInt(usize)) {
        len += 1;
        out_contact_point1.*[1] = contact_points1[point3_index];
        out_contact_point2.*[1] = contact_points2[point3_index];
    }

    out_contact_point1.*[2] = contact_points1[point2_index];
    out_contact_point2.*[2] = contact_points2[point2_index];
    if (point4_index != std.math.maxInt(usize)) {
        len += 1;
        out_contact_point1.*[3] = contact_points1[point4_index];
        out_contact_point2.*[3] = contact_points2[point4_index];
    }
    out_contact_point1.*.len = len;
    out_contact_point2.*.len = len;
}
