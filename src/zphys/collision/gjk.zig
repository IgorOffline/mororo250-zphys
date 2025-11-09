const std = @import("std");
const math = @import("math");

// https://box2d.org/files/ErinCatto_GJK_GDC2010.pdf -> Erin cato presentation for Gjk get closest distance explanation
pub const GjkBox = struct {
    center: math.Vec3,
    orientation: math.Quat,
    half_extents: math.Vec3,

    pub fn support(self: *const @This(), direction: math.Vec3) math.Vec3 {
        if (direction.len2() == 0.0) {
            const local_support = math.vec3(self.half_extents.x(), self.half_extents.y(), self.half_extents.z());
            const world_support = local_support.mulQuat(&self.orientation);
            return self.center.add(&world_support);
        }

        // Localize the direction using the inverse (conjugate) rotation
        const inv_rot = self.orientation.conjugate();
        const local_dir = direction.mulQuat(&inv_rot);

        const sx = if (local_dir.x() >= 0.0) self.half_extents.x() else -self.half_extents.x();
        const sy = if (local_dir.y() >= 0.0) self.half_extents.y() else -self.half_extents.y();
        const sz = if (local_dir.z() >= 0.0) self.half_extents.z() else -self.half_extents.z();

        const local_support = math.vec3(sx, sy, sz);
        const world_support = local_support.mulQuat(&self.orientation);
        return self.center.add(&world_support);
    }

    // Returns the 4 vertices of the face most aligned with the given direction
    // Vertices are returned in counter-clockwise order when viewed from outside
    pub fn getSupportFace(self: *const @This(), direction: math.Vec3) [4]math.Vec3 {
        var face: [4]math.Vec3 = undefined;

        // Transform direction to local space to work with AABB
        const inv_rot = self.orientation.conjugate();
        const local_dir = direction.mulQuat(&inv_rot);

        const abs_x = @abs(local_dir.x());
        const abs_y = @abs(local_dir.y());
        const abs_z = @abs(local_dir.z());

        const hx = self.half_extents.x();
        const hy = self.half_extents.y();
        const hz = self.half_extents.z();

        if (abs_x > abs_y and abs_x > abs_z) {
            const x = if (local_dir.x() > 0) hx else -hx;
            face[0] = math.vec3(x, -hy, -hz);
            face[1] = math.vec3(x, -hy, hz);
            face[2] = math.vec3(x, hy, hz);
            face[3] = math.vec3(x, hy, -hz);
        } else if (abs_y > abs_z) {
            const y = if (local_dir.y() > 0) hy else -hy;
            face[0] = math.vec3(-hx, y, -hz);
            face[1] = math.vec3(-hx, y, hz);
            face[2] = math.vec3(hx, y, hz);
            face[3] = math.vec3(hx, y, -hz);
        } else {
            const z = if (local_dir.z() > 0) hz else -hz;
            face[0] = math.vec3(-hx, -hy, z);
            face[1] = math.vec3(-hx, hy, z);
            face[2] = math.vec3(hx, hy, z);
            face[3] = math.vec3(hx, -hy, z);
        }

        // Transform all vertices from local space to world space
        inline for (0..4) |i| {
            const vertex = face[i].mulQuat(&self.orientation);
            face[i] = self.center.add(&vertex);
        }

        return face;
    }
};

// Algorithm implementation: https://www.youtube.com/watch?v=ajv46BSqcK4&t=887s￥
pub fn gjkIntersect(
    simplex_arrays: [3][]math.Vec3, // [0]=minkowski simplex (A-B), [1]=shape A support points, [2]=shape B support points
    shape_a: anytype,
    shape_b: anytype,
) bool {
    var simplex_size: usize = 0;
    // Keep original naming in code: alias the arrays to descriptive locals
    var simplex = simplex_arrays[0];
    var shape_a_points = simplex_arrays[1];
    var shape_b_points = simplex_arrays[2];

    var search_direction = shape_b.center.sub(&shape_a.center);
    if (search_direction.len2() < 1e-8)
        return true;

    // Add first point from Minkowski difference: S(d) = supportA(d) - supportB(-d)
    const support_a0 = shape_a.support(search_direction);
    const support_b0 = shape_b.support(search_direction.negate());
    simplex[0] = support_a0.sub(&support_b0);
    shape_a_points[0] = support_a0;
    shape_b_points[0] = support_b0;
    simplex_size = 1;

    if (simplex[0].dot(&search_direction) <= 0)
        return false;
    search_direction = simplex[0].negate();

    var iteration: usize = 0;
    while (iteration < 30) : (iteration += 1) {
        const support_a = shape_a.support(search_direction);
        const support_b = shape_b.support(search_direction.negate());
        const new_point = support_a.sub(&support_b);
        if (new_point.dot(&search_direction) <= 0) return false;

        simplex[simplex_size] = new_point;
        shape_a_points[simplex_size] = support_a;
        shape_b_points[simplex_size] = support_b;
        simplex_size += 1;

        const contains_origin = handleSimplex(simplex, shape_a_points, shape_b_points, &simplex_size, &search_direction);
        if (contains_origin) return true;
    }
    return false;
}

fn handleSimplex(simplex: []math.Vec3, shape_a_points: []math.Vec3, shape_b_points: []math.Vec3, simplex_size: *usize, search_direction: *math.Vec3) bool {
    // Based on GJK in 3D - cases line, triangle, tetrahedron
    switch (simplex_size.*) {
        2 => {
            // Line AB (A = last)
            const last_point = simplex[1];
            const previous_point = simplex[0];
            const to_origin = last_point.negate();
            const ab_edge = previous_point.sub(&last_point);

            // New direction perpendicular to AB towards origin
            const ab_cross_ao = ab_edge.cross(&to_origin);
            search_direction.* = ab_cross_ao.cross(&ab_edge);
            if (search_direction.len2() < 1e-12) {
                // pick any perpendicular
                search_direction.* = math.vec3(-ab_edge.y(), ab_edge.x(), 0);
            }
            return false;
        },
        3 => {
            // Triangle ABC (A = last)
            const last_point = simplex[2];
            const point_b = simplex[1];
            const point_c = simplex[0];

            // Capture corresponding A/B support points for reordering
            const support_a_A = shape_a_points[2];
            const support_a_B = shape_a_points[1];
            const support_a_C = shape_a_points[0];
            const support_b_A = shape_b_points[2];
            const support_b_B = shape_b_points[1];
            const support_b_C = shape_b_points[0];

            const to_origin = last_point.negate();
            const ab_edge = point_b.sub(&last_point);
            const ac_edge = point_c.sub(&last_point);
            const triangle_normal = ab_edge.cross(&ac_edge);

            // Determine which side of triangle the origin lies
            const ab_perp_direction = triangle_normal.cross(&ac_edge);
            if (ab_perp_direction.dot(&to_origin) > 0) {
                // Origin is outside AC edge
                simplex[0] = point_c;
                simplex[1] = last_point;
                shape_a_points[0] = support_a_C;
                shape_a_points[1] = support_a_A;
                shape_b_points[0] = support_b_C;
                shape_b_points[1] = support_b_A;

                simplex_size.* = 2;
                search_direction.* = ac_edge.cross(&to_origin).cross(&ac_edge);
                if (search_direction.len2() < 1e-12) search_direction.* = math.vec3(-ac_edge.y(), ac_edge.x(), 0);
                return false;
            }

            const ac_perp_direction = ab_edge.cross(&triangle_normal);
            if (ac_perp_direction.dot(&to_origin) > 0) {
                // Outside AB edge
                simplex[0] = point_b;
                simplex[1] = last_point;
                shape_a_points[0] = support_a_B;
                shape_a_points[1] = support_a_A;
                shape_b_points[0] = support_b_B;
                shape_b_points[1] = support_b_A;

                simplex_size.* = 2;
                search_direction.* = ab_edge.cross(&to_origin).cross(&ab_edge);
                if (search_direction.len2() < 1e-12) search_direction.* = math.vec3(-ab_edge.y(), ab_edge.x(), 0);
                return false;
            }

            // Otherwise, origin is above/below triangle
            if (triangle_normal.dot(&to_origin) > 0) {
                search_direction.* = triangle_normal;
            } else {
                // Wind triangle the other way
                simplex[0] = point_b;
                simplex[1] = point_c;
                simplex[2] = last_point;

                shape_a_points[0] = support_a_B;
                shape_a_points[1] = support_a_C;
                shape_a_points[2] = support_a_A;
                shape_b_points[0] = support_b_B;
                shape_b_points[1] = support_b_C;
                shape_b_points[2] = support_b_A;

                search_direction.* = triangle_normal.negate();
            }
            return false;
        },
        4 => {
            // Tetrahedron ABCD (A = last)
            const last_point = simplex[3];
            const point_b = simplex[2];
            const point_c = simplex[1];
            const point_d = simplex[0];

            // Capture corresponding A/B support points for reordering
            const support_a_A = shape_a_points[3];
            const support_a_B = shape_a_points[2];
            const support_a_C = shape_a_points[1];
            const support_a_D = shape_a_points[0];
            const support_b_A = shape_b_points[3];
            const support_b_B = shape_b_points[2];
            const support_b_C = shape_b_points[1];
            const support_b_D = shape_b_points[0];

            const to_origin = last_point.negate();
            const ab_edge = point_b.sub(&last_point);
            const ac_edge = point_c.sub(&last_point);
            const ad_edge = point_d.sub(&last_point);

            const face_abc = ab_edge.cross(&ac_edge);
            const face_acd = ac_edge.cross(&ad_edge);
            const face_adb = ad_edge.cross(&ab_edge);

            if (face_abc.dot(&to_origin) > 0) {
                simplex[0] = point_c;
                simplex[1] = point_b;
                simplex[2] = last_point;

                shape_a_points[0] = support_a_C;
                shape_a_points[1] = support_a_B;
                shape_a_points[2] = support_a_A;
                shape_b_points[0] = support_b_C;
                shape_b_points[1] = support_b_B;
                shape_b_points[2] = support_b_A;

                simplex_size.* = 3;
                search_direction.* = face_abc;
                return false;
            }
            if (face_acd.dot(&to_origin) > 0) {
                simplex[0] = point_d;
                simplex[1] = point_c;
                simplex[2] = last_point;

                shape_a_points[0] = support_a_D;
                shape_a_points[1] = support_a_C;
                shape_a_points[2] = support_a_A;
                shape_b_points[0] = support_b_D;
                shape_b_points[1] = support_b_C;
                shape_b_points[2] = support_b_A;

                simplex_size.* = 3;
                search_direction.* = face_acd;
                return false;
            }
            if (face_adb.dot(&to_origin) > 0) {
                simplex[0] = point_b;
                simplex[1] = point_d;
                simplex[2] = last_point;

                shape_a_points[0] = support_a_B;
                shape_a_points[1] = support_a_D;
                shape_a_points[2] = support_a_A;
                shape_b_points[0] = support_b_B;
                shape_b_points[1] = support_b_D;
                shape_b_points[2] = support_b_A;

                simplex_size.* = 3;
                search_direction.* = face_adb;
                return false;
            }
            return true;
        },
        else => return false,
    }
}

inline fn signf(value: f32) f32 {
    return if (value >= 0) 1.0 else -1.0;
}

test "gjkIntersect.box_box.separated" {
    const a = GjkBox{
        .center = math.vec3(0, 0, 0),
        .orientation = math.Quat.identity(),
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    const b = GjkBox{
        .center = math.vec3(2.0, 0, 0),
        .orientation = math.Quat.identity(),
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    var simplex_points: [4]math.Vec3 = undefined;
    var shape_a_points: [4]math.Vec3 = undefined;
    var shape_b_points: [4]math.Vec3 = undefined;
    const simplex_arrays: [3][]math.Vec3 = .{ simplex_points[0..], shape_a_points[0..], shape_b_points[0..] };
    try std.testing.expect(!gjkIntersect(simplex_arrays, a, b));
    try std.testing.expect(!gjkIntersect(simplex_arrays, b, a));
}

test "gjkIntersect.box_box.overlap" {
    const a = GjkBox{
        .center = math.vec3(0, 0, 0),
        .orientation = math.Quat.identity(),
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    const b = GjkBox{
        .center = math.vec3(0.75, 0, 0),
        .orientation = math.Quat.identity(),
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    var simplex_points: [4]math.Vec3 = undefined;
    var shape_a_points: [4]math.Vec3 = undefined;
    var shape_b_points: [4]math.Vec3 = undefined;
    const simplex_arrays: [3][]math.Vec3 = .{ simplex_points[0..], shape_a_points[0..], shape_b_points[0..] };
    try std.testing.expect(gjkIntersect(simplex_arrays, a, b));
    try std.testing.expect(gjkIntersect(simplex_arrays, b, a));
}

test "gjkIntersect.box_box.separated_rot_y" {
    const a = GjkBox{
        .center = math.vec3(0, 0, 0),
        .orientation = math.Quat.identity(),
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    const rot = math.Quat.fromAxisAngle(math.vec3(0, 1, 0), math.degreesToRadians(45.0));
    const b = GjkBox{
        .center = math.vec3(2.0, 0, 0),
        .orientation = rot,
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    var simplex_points: [4]math.Vec3 = undefined;
    var shape_a_points: [4]math.Vec3 = undefined;
    var shape_b_points: [4]math.Vec3 = undefined;
    const simplex_arrays: [3][]math.Vec3 = .{ simplex_points[0..], shape_a_points[0..], shape_b_points[0..] };
    try std.testing.expect(!gjkIntersect(simplex_arrays, a, b));
    try std.testing.expect(!gjkIntersect(simplex_arrays, b, a));
}

test "gjkIntersect.box_box.overlap_rot_y" {
    const a = GjkBox{
        .center = math.vec3(0, 0, 0),
        .orientation = math.Quat.identity(),
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    const rot = math.Quat.fromAxisAngle(math.vec3(0, 1, 0), math.degreesToRadians(30.0));
    const b = GjkBox{
        .center = math.vec3(0.5, 0, 0),
        .orientation = rot,
        .half_extents = math.vec3(0.5, 0.5, 0.5),
    };
    var simplex_points: [4]math.Vec3 = undefined;
    var shape_a_points: [4]math.Vec3 = undefined;
    var shape_b_points: [4]math.Vec3 = undefined;
    const simplex_arrays: [3][]math.Vec3 = .{ simplex_points[0..], shape_a_points[0..], shape_b_points[0..] };
    try std.testing.expect(gjkIntersect(simplex_arrays, a, b));
    try std.testing.expect(gjkIntersect(simplex_arrays, b, a));
}
