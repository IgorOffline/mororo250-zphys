const std = @import("std");
const math = @import("math");
const bary = @import("barycentric.zig");

// Vectors have an 16bytes alignment
// would I get any performance different by passing the vectors as reference instead?
pub const EpaResult = struct {
    penetration_axis: math.Vec3,
    penetration_depth: f32,
    collision_point_a: math.Vec3,
    collision_point_b: math.Vec3,
};

// EPA capacity notes:
// - The maximal number of faces for convex shape is F ≤ 2V − 4.
//   Derivation (Euler + edge/face relation):
//   • Euler for convex polyhedra: V − E + F = 2
//   • If all faces are triangles and each edge is shared by exactly two faces: 3F = 2E ⇒ E = 3F/2
//   • Substitute into Euler: V − (3F/2) + F = 2 ⇒ V − F/2 = 2 ⇒ F = 2V − 4
// - For box–box, the CSO (Minkowski difference) can have up to V = 16 vertices,
//   therefore worst-case faces is F ≤ 2*16 − 4 = 28.
// - We use 28 here to safely handle box–box with a Minkowski buffer of 16 vertices.
//   Consider deriving this from the supplied buffer length in the future.
const simplex_size = 4;
const tolerance = 0.001;
const max_faces = 28;

const Edge = struct { a: u32, b: u32 };

/// Algorithm reference: https://winter.dev/articles/epa-algorithm
/// Accepts array-of-arrays matching GJK layout:
///  - simplex_arrays[0]: Minkowski points (A - B)
///  - simplex_arrays[1]: shape A support points
///  - simplex_arrays[2]: shape B support points
pub fn epa(simplex_arrays: [3][]math.Vec3, shape_a: anytype, shape_b: anytype) EpaResult {
    var polytype = simplex_arrays[0];
    var shape_a_points = simplex_arrays[1];
    var shape_b_points = simplex_arrays[2];

    const max_vertices: u32 = @intCast(polytype.len);
    // Guard: for a triangulated convex polytope with V vertices, the maximal number of faces is F <= 2V - 4
    std.debug.assert(max_faces >= (2 * max_vertices - 4));

    // Face index buffer large enough for worst-case faces for box - box (28)
    var face_edge_indexes: [max_faces * 3]u32 = undefined;
    face_edge_indexes[0..12].* = [12]u32{ 0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2 };
    var normals: [max_faces]math.Vec3 = undefined;
    var distances: [max_faces]f32 = undefined;
    var centroids: [max_faces]math.Vec3 = undefined;
    var horizons: [max_faces]Edge = undefined;
    var min_index: u32 = 0;
    var min_distance: f32 = std.math.floatMax(f32);

    comptime var i: u32 = 0;
    inline while (i < simplex_size) : (i += 1) {
        const face = face_edge_indexes[i * 3 ..][0..3];
        const c = triangleCentroid(polytype, face);
        centroids[i] = c;
        calcNormalDistance(polytype, face, c, &normals[i], &distances[i]);

        if (distances[i] < min_distance) {
            min_distance = distances[i];
            min_index = i;
        }
    }

    const max_number_iter = 30;
    var iter: u32 = 0;
    var edges_count: u32 = simplex_size;
    var face_count: u32 = edges_count;
    while (iter < max_number_iter) : (iter += 1) {
        const min_normal = normals[min_index];
        min_distance = distances[min_index];

        const support_a = shape_a.support(min_normal);
        const support_b = shape_b.support(min_normal.negate());
        const support = support_a.sub(&support_b);
        const support_distance: f32 = min_normal.dot(&support);

        // Convergence
        if (support_distance - min_distance < tolerance) {
            const face = face_edge_indexes[min_index * 3 ..][0..3];
            const a_idx = face[0];
            const b_idx = face[1];
            const c_idx = face[2];

            const closest_point = min_normal.mulScalar(min_distance);
            const weights = bary.barycentricTriangle(closest_point, polytype[a_idx], polytype[b_idx], polytype[c_idx]);

            const collision_point_a =
                shape_a_points[a_idx].mulScalar(weights.x())
                    .add(&shape_a_points[b_idx].mulScalar(weights.y()))
                    .add(&shape_a_points[c_idx].mulScalar(weights.z()));
            const collision_point_b =
                shape_b_points[a_idx].mulScalar(weights.x())
                    .add(&shape_b_points[b_idx].mulScalar(weights.y()))
                    .add(&shape_b_points[c_idx].mulScalar(weights.z()));

            return .{
                .penetration_axis = min_normal.mulScalar(min_distance),
                .penetration_depth = min_distance,
                .collision_point_a = collision_point_a,
                .collision_point_b = collision_point_b,
            };
        }

        // Expand polytope: remove visible faces, collect horizon edges
        var j: u32 = 0;
        var new_face_count: u32 = 0;
        var horizon_count: u32 = 0;
        min_distance = std.math.floatMax(f32);
        while (j < face_count) : (j += 1) {
            const face = face_edge_indexes[j * 3 ..][0..3];
            const centroid = centroids[j];
            // Using centroid here improves numerical stability
            const diff = support.sub(&centroid);

            if (normals[j].dot(&diff) > 0) {
                buildHorizon(horizons[0..], &horizon_count, face);
                continue;
            }

            // Keep this face: compact into new_face_count slot
            if (j != new_face_count) {
                var dst: *[3]u32 = face_edge_indexes[new_face_count * 3 ..][0..3];
                dst[0] = face[0];
                dst[1] = face[1];
                dst[2] = face[2];
            }
            normals[new_face_count] = normals[j];
            distances[new_face_count] = distances[j];
            centroids[new_face_count] = centroids[j];
            if (distances[new_face_count] < min_distance) {
                min_distance = distances[new_face_count];
                min_index = new_face_count;
            }

            new_face_count += 1;
        }
        face_count = new_face_count;

        // Build new faces using the horizon
        reconstructFaces(face_edge_indexes[0..], &face_count, horizons[0..], horizon_count, edges_count);
        std.debug.assert(edges_count < polytype.len);
        polytype[edges_count] = support;
        shape_a_points[edges_count] = support_a;
        shape_b_points[edges_count] = support_b;
        edges_count += 1;

        // Recompute normals/distances for new faces only
        j = new_face_count;
        while (j < face_count) : (j += 1) {
            const face = face_edge_indexes[j * 3 ..][0..3];
            const c = triangleCentroid(polytype, face);
            centroids[j] = c;

            calcNormalDistance(polytype, face, c, &normals[j], &distances[j]);
            if (min_distance > distances[j]) {
                min_distance = distances[j];
                min_index = j;
            }
        }
    }

    const face = face_edge_indexes[min_index * 3 ..][0..3];
    const a_idx = face[0];
    const b_idx = face[1];
    const c_idx = face[2];

    const closest_point = normals[min_index].mulScalar(min_distance);
    const weights = bary.barycentricTriangle(closest_point, polytype[a_idx], polytype[b_idx], polytype[c_idx]);

    const collision_point_a =
        shape_a_points[a_idx].mulScalar(weights.x())
            .add(&shape_a_points[b_idx].mulScalar(weights.y()))
            .add(&shape_a_points[c_idx].mulScalar(weights.z()));
    const collision_point_b =
        shape_b_points[a_idx].mulScalar(weights.x())
            .add(&shape_b_points[b_idx].mulScalar(weights.y()))
            .add(&shape_b_points[c_idx].mulScalar(weights.z()));

    return .{
        .penetration_axis = normals[min_index],
        .penetration_depth = min_distance,
        .collision_point_a = collision_point_a,
        .collision_point_b = collision_point_b,
    };
}

// Inline helpers (no copies, preserve logic/style)
inline fn triangleCentroid(polytype: []math.Vec3, face: []const u32) math.Vec3 {
    const p0 = polytype[face[0]];
    const p1 = polytype[face[1]];
    const p2 = polytype[face[2]];
    return p0.add(&p1).add(&p2).mulScalar(1.0 / 3.0);
}
inline fn calcNormalDistance(polytype: []math.Vec3, face: []const u32, centroid: math.Vec3, normal: *math.Vec3, distance: *f32) void {
    // Use centroid as the local origin to improve numerical stability for far-away triangles
    const a = polytype[face[1]].sub(&centroid);
    const b = polytype[face[2]].sub(&centroid);
    normal.* = a.cross(&b).normalize(math.eps_f32);
    distance.* = normal.*.dot(&centroid);

    // Ensure normal points away from the origin (0,0,0)
    if (distance.* < 0) {
        distance.* *= -1;
        normal.* = normal.*.negate();
    }
}

inline fn buildHorizon(buffer: []Edge, len: *u32, face: []const u32) void {
    // Face is visible from the new support point -> collect its edges for the horizon
    addIfUnique(buffer, len, .{ .a = face[0], .b = face[1] });
    addIfUnique(buffer, len, .{ .a = face[1], .b = face[2] });
    addIfUnique(buffer, len, .{ .a = face[2], .b = face[0] });
}

inline fn addIfUnique(buffer: []Edge, len: *u32, value: Edge) void {
    // Single-pass: cancel reverse, skip exact duplicate, else insert
    var i: u32 = 0;
    while (i < len.*) : (i += 1) {
        const e = buffer[i];
        // Reverse edge found -> cancel (symmetric difference)
        if (e.a == value.b and e.b == value.a) {
            len.* -= 1;
            buffer[i] = buffer[len.*]; // swap-erase
            return;
        }
        // Exact same oriented edge -> duplicate, skip
        if (e.a == value.a and e.b == value.b) {
            return;
        }
    }
    // Insert new boundary edge
    buffer[len.*] = value;
    len.* += 1;
}

inline fn reconstructFaces(face_edge_indexes: []u32, face_count: *u32, horizons: []const Edge, horizon_count: u32, edges_count: u32) void {
    // Build new faces using the horizon
    std.debug.assert(face_count.* + horizon_count <= max_faces);
    for (horizons[0..horizon_count]) |horizon| {
        const new_face = face_edge_indexes[face_count.* * 3 ..][0..3];
        new_face[0] = horizon.a;
        new_face[1] = horizon.b;
        new_face[2] = edges_count;
        face_count.* += 1;
    }
}
