const std = @import("std");
const math = @import("math");

pub const EpaResult = struct {
    normal: math.Vec3,
    penetration_depth: f32,
};

// For now I multiply the size for 2 of the arrays to support Minkowski difference of 2 quads which may have max size of 8 vertices and 12 faces
// we could get this value as a multiplication of two comptime sizes fromm shape_a and shape_b
const simplex_size = 4;
const tolerance = 0.001;
const max_faces = 12;

const Edge = struct { a: u32, b: u32 };

// Algorithm reference: https://winter.dev/articles/epa-algorithm
pub fn epa(polytype: [] math.Vec3, shape_a: anytype, shape_b: anytype) EpaResult {
    // Face index buffer large enough for worst-case faces (12)
    var face_edge_indexes: [max_faces * 3]u32 = undefined;
    face_edge_indexes[0..12].* = [12]u32{ 0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2 };
    var normals: [max_faces]math.Vec3 = undefined;
    var distances: [max_faces]f32 = undefined;
    var horizons: [max_faces]Edge = undefined;
    var min_index: u32 = 0;
    var min_distance: f32 = std.math.floatMax(f32);

    comptime var i: u32 = 0;
    inline while (i < simplex_size) : (i += 1) {
        const face = face_edge_indexes[i * 3 ..][0..3];
        calcNormalDistance(polytype, face, &normals[i], &distances[i]);

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

        const support = shape_a.support(min_normal).sub(&shape_b.support(min_normal.negate()));
        const support_distance: f32 = min_normal.dot(&support);

        // Convergence
        if (support_distance - min_distance < tolerance) {
            return .{ .normal = min_normal, .penetration_depth = min_distance };
        }

        // Expand polytope: remove visible faces, collect horizon edges
        var j: u32 = 0;
        var new_face_count: u32 = 0;
        var horizon_count: u32 = 0;
        min_distance = std.math.floatMax(f32);
        while (j < face_count) : (j += 1) {
            const face = face_edge_indexes[j * 3 ..][0..3];
            const diff = support.sub(&polytype[face[0]]);

            if (normals[j].dot(&diff) > 0) {
                buildHorizon(horizons[0..], &horizon_count, face);
                continue;
            }

            if (j == new_face_count) continue;
            var dst : *[3]u32 = face_edge_indexes[face_count * 3..][0..3];
            dst[0] = face[0];
            dst[1] = face[1];
            dst[2] = face[2];
            normals[new_face_count] = normals[j];
            distances[new_face_count] = distances[j];
            if (distances[j] < min_distance) {
                min_distance = distances[j];
                min_index = new_face_count;
            }

            new_face_count += 1;
        }
        face_count = new_face_count;

        // Build new faces using the horizon
        reconstructFaces(face_edge_indexes[0..], &face_count, horizons[0..], horizon_count, edges_count);
        std.debug.assert(edges_count < polytype.len);
        polytype[edges_count] = support;
        edges_count += 1;

        // Recompute normals/distances for new faces only
        j = new_face_count;
        while (j < face_count) : (j += 1) {
            const face = face_edge_indexes[j * 3 ..][0..3];

            calcNormalDistance(polytype, face, &normals[j], &distances[j]);
            if (min_distance > distances[j]) {
                min_distance = distances[j];
                min_index = j;
            }
        }
    }

    return .{ .normal = normals[min_index], .penetration_depth = min_distance };
}

// Inline helpers (no copies, preserve logic/style)
inline fn calcNormalDistance(polytype: [] math.Vec3, face: []const u32, normal: *math.Vec3, distance: *f32) void {
    const a = polytype[face[1]].sub(&polytype[face[0]]);
    const b = polytype[face[2]].sub(&polytype[face[0]]);
    normal.* = a.cross(&b).normalize(math.eps_f32);
    distance.* = normal.*.dot(&polytype[face[0]]);

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
    for (buffer[0..len.*]) |item| {
        if (item.a == value.a and item.b == value.b) {
            return;
        }
    }
    buffer[len.*] = value;
    len.* += 1;
}

inline fn reconstructFaces(face_edge_indexes: []u32, face_count: *u32, horizons: []const Edge, horizon_count: u32, edges_count: u32) void {
    // Build new faces using the horizon
    for (horizons[0..horizon_count]) |horizon| {
        const new_face = face_edge_indexes[face_count.* * 3 ..][0..3];
        new_face[0] = horizon.a;
        new_face[1] = horizon.b;
        new_face[2] = edges_count;
        face_count.* += 1;
    }
}
