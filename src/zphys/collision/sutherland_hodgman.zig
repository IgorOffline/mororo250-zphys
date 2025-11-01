const math = @import("math");

/// This function will only work for convex polygon -> for non convex we would need to do an specialized function
/// that returns bigger array
/// Implements sutherland-hodgman algorithm to clip one polygon against another by
/// iteratively clipping against each edge
/// output size n + m -> maximum output size of clipped convex polygon
///     - Each of the clipping planes can, at most, add one new vertex to the convex polygon,
///     - and this process repeats M times, leading to a maximum of N+M vertices.
pub fn clipPolyPoly(comptime max_length: usize,
    poly_to_clip: *const []math.Vec3,
    clipping_poly: *const []math.Vec3,
    clipping_normal : math.Vec3,
    out_poly:*[max_length]math.Vec3) []math.Vec3 {
    var buffer : [max_length]math.Vec3 = undefined;

    // Initialize first buffer with input polygon
    for (poly_to_clip, 0..) |vertex, idx| {
        out_poly[idx] = vertex;
    }
    var current_len: usize = poly_to_clip.len;

    var i:u32 = 0;
    while (i < clipping_poly.len) : (i += 1){
        const vertex1 = clipping_poly[i];
        const vertex2 = clipping_poly[(i + 1) % clipping_poly.len];
        // inward clip normal
        const clip_normal = clipping_normal.cross(vertex2.sub(&vertex1));

        // Double buffering: swap between out_poly and buffer
        if (i % 2 == 0) {
            current_len = clipPolyPlane(out_poly[0..current_len], vertex1, clip_normal, &buffer);
        } else {
            current_len = clipPolyPlane(buffer[0..current_len], vertex1, clip_normal, out_poly);
        }

        if (current_len == 0) {
            return out_poly[0..0];
        }
    }

    if (clipping_poly.len % 2 == 0) {
        for (buffer[0..current_len], 0..) |vertex, idx| {
            out_poly[idx] = vertex;
        }
    }

    return out_poly[0..current_len];
}

fn clipPolyPlane(poly_to_clip: []const math.Vec3, plane_origin: math.Vec3, plane_normal : math.Vec3, out_poly: []math.Vec3) usize {
    var prev_vertex = poly_to_clip[poly_to_clip.len - 1];
    var prev_num = (plane_origin - prev_vertex).dot(plane_normal);
    var prev_inside = prev_num < 0;

    var i:u32 = 0;
    var out_len: u32 = 0;
    while (i < poly_to_clip.len) : (i += 1) {
        const cur_vertex =  poly_to_clip[i];
        const cur_num = (plane_origin - cur_vertex).dot(plane_normal);
        const cur_inside = cur_num < 0;

        if (cur_inside != prev_inside) {
            const cur_prev = cur_vertex - prev_vertex;
            const denom = cur_prev.dot(plane_normal);
            if (denom != 0) {
                out_poly[out_len] = prev_vertex.add(&cur_prev.mulScalar(prev_num/denom));
                out_len += 1;
            }
            else {
                cur_inside = prev_inside; // edge is parallel to plane, treat point as if it were on the same side as the last point
            }
        }

        if (cur_inside){
            out_poly[out_len] = cur_vertex;
            out_len += 1;
            prev_inside = true;
        }

        prev_vertex = cur_vertex;
        prev_num = cur_num;
        prev_inside = cur_inside;
    }
    
    return out_len;
}

const std = @import("std");

test "clipPolyPoly - partial overlap" {
    // Two squares partially overlapping
    const square1 = [_]math.Vec3{
        math.Vec3.new(-1, -1, 0),
        math.Vec3.new(1, -1, 0),
        math.Vec3.new(1, 1, 0),
        math.Vec3.new(-1, 1, 0),
    };
    
    const square2 = [_]math.Vec3{
        math.Vec3.new(0, -1, 0),
        math.Vec3.new(2, -1, 0),
        math.Vec3.new(2, 1, 0),
        math.Vec3.new(0, 1, 0),
    };
    
    const normal = math.Vec3.new(0, 0, 1);
    var out_poly: [8]math.Vec3 = undefined;
    
    const result = clipPolyPoly(4, &square1, 4, &square2, normal, &out_poly);
    
    // The intersection should be a rectangle from (0,-1) to (1,1)
    try std.testing.expect(result.len == 4);
    
    // Check that all resulting vertices are within both polygons
    for (result) |vertex| {
        try std.testing.expect(vertex.x() >= 0 and vertex.x() <= 1);
        try std.testing.expect(vertex.y() >= -1 and vertex.y() <= 1);
    }
}

test "clipPolyPoly - complete containment" {
    // Small square inside larger square
    const large_square = [_]math.Vec3{
        math.Vec3.new(-2, -2, 0),
        math.Vec3.new(2, -2, 0),
        math.Vec3.new(2, 2, 0),
        math.Vec3.new(-2, 2, 0),
    };
    
    const small_square = [_]math.Vec3{
        math.Vec3.new(-1, -1, 0),
        math.Vec3.new(1, -1, 0),
        math.Vec3.new(1, 1, 0),
        math.Vec3.new(-1, 1, 0),
    };
    
    const normal = math.Vec3.new(0, 0, 1);
    var out_poly: [8]math.Vec3 = undefined;
    
    const result = clipPolyPoly(4, &small_square, 4, &large_square, normal, &out_poly);
    
    // Small square should be completely preserved
    try std.testing.expect(result.len == 4);
    
    // Verify vertices match the small square (may be in different order)
    for (result) |vertex| {
        try std.testing.expect(vertex.x() >= -1 and vertex.x() <= 1);
        try std.testing.expect(vertex.y() >= -1 and vertex.y() <= 1);
    }
}

test "clipPolyPoly - no overlap" {
    // Two squares that don't intersect
    const square1 = [_]math.Vec3{
        math.Vec3.new(-2, -2, 0),
        math.Vec3.new(-1, -2, 0),
        math.Vec3.new(-1, -1, 0),
        math.Vec3.new(-2, -1, 0),
    };
    
    const square2 = [_]math.Vec3{
        math.Vec3.new(1, 1, 0),
        math.Vec3.new(2, 1, 0),
        math.Vec3.new(2, 2, 0),
        math.Vec3.new(1, 2, 0),
    };
    
    const normal = math.Vec3.new(0, 0, 1);
    var out_poly: [8]math.Vec3 = undefined;
    
    const result = clipPolyPoly(4, &square1, 4, &square2, normal, &out_poly);
    
    // No intersection, result should be empty
    try std.testing.expect(result.len == 0);
}

test "clipPolyPoly - triangle clips square" {
    // Square
    const square = [_]math.Vec3{
        math.Vec3.new(-1, -1, 0),
        math.Vec3.new(1, -1, 0),
        math.Vec3.new(1, 1, 0),
        math.Vec3.new(-1, 1, 0),
    };
    
    // Triangle covering top-right portion
    const triangle = [_]math.Vec3{
        math.Vec3.new(0, 0, 0),
        math.Vec3.new(2, 0, 0),
        math.Vec3.new(0, 2, 0),
    };
    
    const normal = math.Vec3.new(0, 0, 1);
    var out_poly: [7]math.Vec3 = undefined;
    
    const result = clipPolyPoly(4, &square, 3, &triangle, normal, &out_poly);
    
    // Should have intersection
    try std.testing.expect(result.len > 0);
    
    // All vertices should be within valid bounds
    for (result) |vertex| {
        try std.testing.expect(vertex.x() >= -1 and vertex.x() <= 2);
        try std.testing.expect(vertex.y() >= -1 and vertex.y() <= 2);
        try std.testing.expect(vertex.z() == 0);
    }
}

test "clipPolyPoly - identical polygons" {
    // Two identical squares
    const square = [_]math.Vec3{
        math.Vec3.new(-1, -1, 0),
        math.Vec3.new(1, -1, 0),
        math.Vec3.new(1, 1, 0),
        math.Vec3.new(-1, 1, 0),
    };
    
    const normal = math.Vec3.new(0, 0, 1);
    var out_poly: [8]math.Vec3 = undefined;
    
    const result = clipPolyPoly(4, &square, 4, &square, normal, &out_poly);
    
    // Result should be the same square
    try std.testing.expect(result.len == 4);
    
    for (result) |vertex| {
        try std.testing.expect(vertex.x() >= -1 and vertex.x() <= 1);
        try std.testing.expect(vertex.y() >= -1 and vertex.y() <= 1);
    }
}
