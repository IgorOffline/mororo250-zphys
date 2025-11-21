const std = @import("std");
const math = @import("math");
const zphys = @import("zphys");
const rl = @import("raylib");
const DebugRenderer = @import("shared/debug_renderer.zig").DebugRenderer;
const SceneRenderer = @import("shared/scene_renderer.zig").SceneRenderer;

// Todo: unify logis of the state and the internal epa and gjk logic.
// Is it possible to do it while keeping pefomance?
const GjkBox = zphys.collision.gjk.GjkBox;

const GJKStepper = struct {
    simplex: [4]math.Vec3 = undefined,
    shape_a_points: [4]math.Vec3 = undefined,
    shape_b_points: [4]math.Vec3 = undefined,
    simplex_size: usize = 0,
    search_direction: math.Vec3 = math.vec3(0, 0, 0),
    iteration: usize = 0,
    finished: bool = false,
    intersect: bool = false,
    
    // Output buffers for EPA
    out_simplex: [16]math.Vec3 = undefined,
    out_shape_a: [16]math.Vec3 = undefined,
    out_shape_b: [16]math.Vec3 = undefined,
    out_count: usize = 0,

    pub fn init(shape_a: GjkBox, shape_b: GjkBox) GJKStepper {
        var self = GJKStepper{};
        self.search_direction = shape_b.center.sub(&shape_a.center);
        
        if (self.search_direction.len2() < 1e-8) {
            self.finished = true;
            self.intersect = true;
        }
        return self;
    }

    pub fn step(self: *GJKStepper, shape_a: GjkBox, shape_b: GjkBox) void {
        if (self.finished or self.iteration >= 32) {
            self.finished = true;
            return;
        }
        self.iteration += 1;

        const result = zphys.collision.gjk.gjkStep(
            self.simplex[0..],
            self.shape_a_points[0..],
            self.shape_b_points[0..],
            &self.simplex_size,
            &self.search_direction,
            shape_a,
            shape_b,
        );

        switch (result) {
            .Intersection => {
                self.finished = true;
                self.intersect = true;
                // Copy simplex to output buffers for EPA
                for (0..self.simplex_size) |i| {
                    self.out_simplex[i] = self.simplex[i];
                    self.out_shape_a[i] = self.shape_a_points[i];
                    self.out_shape_b[i] = self.shape_b_points[i];
                }
                self.out_count = self.simplex_size;
            },
            .NoIntersection => {
                self.finished = true;
                self.intersect = false;
            },
            .Continue => {},
        }
    }
};

const Edge = struct { a: u32, b: u32 };
const max_faces = 28;
const simplex_size = 4;
const tolerance = 0.001;

const EPAState = struct {
    polytype: [16]math.Vec3 = undefined,
    shape_a_points: [16]math.Vec3 = undefined,
    shape_b_points: [16]math.Vec3 = undefined,
    
    face_edge_indexes: [max_faces * 3]u32 = undefined,
    normals: [max_faces]math.Vec3 = undefined,
    distances: [max_faces]f32 = undefined,
    centroids: [max_faces]math.Vec3 = undefined,
    horizons: [max_faces]Edge = undefined,
    
    // State vars
    min_index: u32 = 0,
    min_distance: f32 = std.math.floatMax(f32),
    iter: u32 = 0,
    edges_count: u32 = 0,
    face_count: u32 = 0,
    
    initialized: bool = false,
    finished: bool = false,
    
    // Debug - current support info
    current_support: math.Vec3 = math.vec3(0,0,0),
    current_min_normal: math.Vec3 = math.vec3(0,0,0),

    pub fn init(gjk_state: GJKStepper) EPAState {
        var self = EPAState{};
        
        // Copy from GJK output
        for (0..gjk_state.out_count) |i| {
            self.polytype[i] = gjk_state.out_simplex[i];
            self.shape_a_points[i] = gjk_state.out_shape_a[i];
            self.shape_b_points[i] = gjk_state.out_shape_b[i];
        }
        self.edges_count = @intCast(gjk_state.out_count);
        self.face_count = self.edges_count;
        
        // Initial Faces
        self.face_edge_indexes[0..12].* = [12]u32{ 0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2 };
        
        // Initial Normals
        var i: u32 = 0;
        while (i < 4) : (i += 1) {
            const face = self.face_edge_indexes[i * 3 ..][0..3];
            const c = triangleCentroid(self.polytype[0..], face);
            self.centroids[i] = c;
            calcNormalDistance(self.polytype[0..], face, c, &self.normals[i], &self.distances[i]);

            if (self.distances[i] < self.min_distance) {
                self.min_distance = self.distances[i];
                self.min_index = i;
            }
        }
        self.initialized = true;
        return self;
    }

    pub fn step(self: *EPAState, shape_a: GjkBox, shape_b: GjkBox) void {
        if (self.finished or self.iter >= 30) {
            self.finished = true;
            std.log.info("EPA: Max iter reached", .{});
            return;
        }
        self.iter += 1;

        const min_normal = self.normals[self.min_index];
        self.min_distance = self.distances[self.min_index];
        self.current_min_normal = min_normal;

        const support_a = shape_a.support(min_normal);
        const support_b = shape_b.support(min_normal.negate());
        const support = support_a.sub(&support_b);
        self.current_support = support;
        
        const support_distance: f32 = min_normal.dot(&support);

        if (support_distance - self.min_distance < tolerance) {
            self.finished = true;
            std.log.info("EPA: Converged! Iter={}", .{self.iter});
            return;
        }

        var j: u32 = 0;
        var new_face_count: u32 = 0;
        var horizon_count: u32 = 0;
        self.min_distance = std.math.floatMax(f32);

        while (j < self.face_count) : (j += 1) {
            const face = self.face_edge_indexes[j * 3 ..][0..3];
            const centroid = self.centroids[j];
            const diff = support.sub(&centroid);

            if (self.normals[j].dot(&diff) > 0) {
                buildHorizon(self.horizons[0..], &horizon_count, face);
                continue;
            }

            if (j != new_face_count) {
                var dst = self.face_edge_indexes[new_face_count * 3 ..][0..3];
                dst[0] = face[0];
                dst[1] = face[1];
                dst[2] = face[2];
            }
            self.normals[new_face_count] = self.normals[j];
            self.distances[new_face_count] = self.distances[j];
            self.centroids[new_face_count] = self.centroids[j];
            if (self.distances[new_face_count] < self.min_distance) {
                self.min_distance = self.distances[new_face_count];
                self.min_index = new_face_count;
            }

            new_face_count += 1;
        }
        self.face_count = new_face_count;

        reconstructFaces(self.face_edge_indexes[0..], &self.face_count, self.horizons[0..], horizon_count, self.edges_count);
        
        if (self.edges_count < self.polytype.len) {
            self.polytype[self.edges_count] = support;
            self.shape_a_points[self.edges_count] = support_a;
            self.shape_b_points[self.edges_count] = support_b;
            self.edges_count += 1;
        } else {
            std.log.warn("EPA: Buffer full! Edges count {}", .{self.edges_count});
            self.finished = true; // Or crash if we want to simulate exact behavior (assertion failure)
        }

        j = new_face_count;
        while (j < self.face_count) : (j += 1) {
            const face = self.face_edge_indexes[j * 3 ..][0..3];
            const c = triangleCentroid(self.polytype[0..], face);
            self.centroids[j] = c;

            calcNormalDistance(self.polytype[0..], face, c, &self.normals[j], &self.distances[j]);
            if (self.min_distance > self.distances[j]) {
                self.min_distance = self.distances[j];
                self.min_index = j;
            }
        }
    }
};

// Helpers copied from epa.zig
fn triangleCentroid(polytype: []math.Vec3, face: []const u32) math.Vec3 {
    const p0 = polytype[face[0]];
    const p1 = polytype[face[1]];
    const p2 = polytype[face[2]];
    return p0.add(&p1).add(&p2).mulScalar(1.0 / 3.0);
}
fn calcNormalDistance(polytype: []math.Vec3, face: []const u32, centroid: math.Vec3, normal: *math.Vec3, distance: *f32) void {
    const a = polytype[face[1]].sub(&centroid);
    const b = polytype[face[2]].sub(&centroid);
    normal.* = a.cross(&b).normalize(math.eps_f32);
    distance.* = normal.*.dot(&centroid);
    if (distance.* < 0) {
        distance.* *= -1;
        normal.* = normal.*.negate();
    }
}
fn buildHorizon(buffer: []Edge, len: *u32, face: []const u32) void {
    addIfUnique(buffer, len, .{ .a = face[0], .b = face[1] });
    addIfUnique(buffer, len, .{ .a = face[1], .b = face[2] });
    addIfUnique(buffer, len, .{ .a = face[2], .b = face[0] });
}
fn addIfUnique(buffer: []Edge, len: *u32, value: Edge) void {
    var i: u32 = 0;
    while (i < len.*) : (i += 1) {
        const e = buffer[i];
        if (e.a == value.b and e.b == value.a) {
            len.* -= 1;
            buffer[i] = buffer[len.*];
            return;
        }
        if (e.a == value.a and e.b == value.b) return;
    }
    buffer[len.*] = value;
    len.* += 1;
}
fn reconstructFaces(face_edge_indexes: []u32, face_count: *u32, horizons: []const Edge, horizon_count: u32, edges_count: u32) void {
    const max_faces_local = 28;
    std.debug.assert(face_count.* + horizon_count <= max_faces_local);
    for (horizons[0..horizon_count]) |horizon| {
        const new_face = face_edge_indexes[face_count.* * 3 ..][0..3];
        new_face[0] = horizon.a;
        new_face[1] = horizon.b;
        new_face[2] = edges_count;
        face_count.* += 1;
    }
}

pub fn main() !void {
    const screenWidth = 800;
    const screenHeight = 450;
    rl.initWindow(screenWidth, screenHeight, "GJK/EPA Visualizer");
    defer rl.closeWindow();

    var camera = rl.Camera{
        .position = .init(5, 5, 5),
        .target = .init(0, 0, 0),
        .up = .init(0, 1, 0),
        .fovy = 45,
        .projection = .perspective,
    };

    rl.setTargetFPS(60);
    rl.disableCursor();

    // Test Case: Ground vs Box (Exact touching) - Replicates Pyramid crash scenario (first pair)
    // Ground: Center (0, -1, 0), HalfExtents (50, 1, 50) (Top at Y=0)
    const box_a = GjkBox{ .center = math.vec3(0, -1, 0), .orientation = math.Quat.identity(), .half_extents = math.vec3(50, 1, 50) };
    // Box: Center (-4.725, 0.5, 0), HalfExtents (0.5, 0.5, 0.5) (Bottom at Y=0) - Matches id=1 from pyramid
    const box_b = GjkBox{ .center = math.vec3(-4.725, 0.5, 0), .orientation = math.Quat.identity(), .half_extents = math.vec3(0.5, 0.5, 0.5) };

    var gjk_state = GJKStepper.init(box_a, box_b);
    var epa_state: ?EPAState = null;
    
    var paused = true;
    var auto_step_timer: f32 = 0;
    const auto_step_delay: f32 = 0.5;

    while (!rl.windowShouldClose()) {
        updateGJK_EPA(&gjk_state, &epa_state, &paused, &auto_step_timer, auto_step_delay, box_a, box_b);
        camera.update(.free);

        rl.beginDrawing();
        defer rl.endDrawing();
        rl.clearBackground(.white);

        rl.beginMode3D(camera);
        drawScene(box_a, box_b, gjk_state, epa_state);
        rl.endMode3D();

        drawUI(gjk_state, epa_state);
    }
}

fn updateGJK_EPA(gjk_state: *GJKStepper, epa_state: *?EPAState, paused: *bool, auto_step_timer: *f32, auto_step_delay: f32, box_a: GjkBox, box_b: GjkBox) void {
    if (rl.isKeyPressed(.r)) {
        gjk_state.* = GJKStepper.init(box_a, box_b);
        epa_state.* = null;
        paused.* = true;
    }
    if (rl.isKeyPressed(.p)) paused.* = !paused.*;

    var do_step = false;
    if (rl.isKeyPressed(.n)) do_step = true;

    if (!paused.*) {
        auto_step_timer.* += rl.getFrameTime();
        if (auto_step_timer.* > auto_step_delay) {
            auto_step_timer.* = 0;
            do_step = true;
        }
    }

    if (do_step) {
        if (!gjk_state.finished) {
             gjk_state.step(box_a, box_b);
             if (gjk_state.finished and gjk_state.intersect) {
                 epa_state.* = EPAState.init(gjk_state.*);
             }
        } else if (epa_state.*) |*epa_s| {
             if (!epa_s.finished) {
                 epa_s.step(box_a, box_b);
             }
        }
    }
}

fn drawScene(box_a: GjkBox, box_b: GjkBox, gjk_state: GJKStepper, epa_state: ?EPAState) void {
    // Draw Origin
    rl.drawGrid(10, 1.0);
    rl.drawLine3D(.init(0,0,0), .init(1,0,0), .red);
    rl.drawLine3D(.init(0,0,0), .init(0,1,0), .green);
    rl.drawLine3D(.init(0,0,0), .init(0,0,1), .blue);

    rl.drawCubeWiresV(.init(box_a.center.x(), box_a.center.y(), box_a.center.z()),
                      .init(box_a.half_extents.x()*2, box_a.half_extents.y()*2, box_a.half_extents.z()*2), .gray);
     rl.drawCubeWiresV(.init(box_b.center.x(), box_b.center.y(), box_b.center.z()),
                      .init(box_b.half_extents.x()*2, box_b.half_extents.y()*2, box_b.half_extents.z()*2), .dark_gray);

    // Draw GJK Simplex
    if (!gjk_state.finished or !gjk_state.intersect) {
        for (0..gjk_state.simplex_size) |i| {
            const p = gjk_state.simplex[i];
            rl.drawSphere(.init(p.x(), p.y(), p.z()), 0.05, .orange);
        }
        // Draw lines between points (naive)
        if (gjk_state.simplex_size > 1) {
            for (0..gjk_state.simplex_size) |i| {
                for (i+1..gjk_state.simplex_size) |j| {
                     const p1 = gjk_state.simplex[i];
                     const p2 = gjk_state.simplex[j];
                     rl.drawLine3D(.init(p1.x(), p1.y(), p1.z()), .init(p2.x(), p2.y(), p2.z()), .orange);
                }
            }
        }
    }

    // Draw EPA Polytope
    if (epa_state) |epa_s| {
        // Draw vertices
        for (0..epa_s.edges_count) |i| {
             const p = epa_s.polytype[i];
             rl.drawSphere(.init(p.x(), p.y(), p.z()), 0.03, .purple);
        }

        // Draw faces (wireframe)
        var f: u32 = 0;
        while (f < epa_s.face_count) : (f += 1) {
             const face = epa_s.face_edge_indexes[f * 3 ..][0..3];
             const p0 = epa_s.polytype[face[0]];
             const p1 = epa_s.polytype[face[1]];
             const p2 = epa_s.polytype[face[2]];
             const v0 = rl.Vector3.init(p0.x(), p0.y(), p0.z());
             const v1 = rl.Vector3.init(p1.x(), p1.y(), p1.z());
             const v2 = rl.Vector3.init(p2.x(), p2.y(), p2.z());

             rl.drawLine3D(v0, v1, .purple);
             rl.drawLine3D(v1, v2, .purple);
             rl.drawLine3D(v2, v0, .purple);

             // Draw normal
             const c = epa_s.centroids[f];
             const n = epa_s.normals[f];
             const vc = rl.Vector3.init(c.x(), c.y(), c.z());
             const ve = rl.Vector3.init(c.x() + n.x()*0.2, c.y() + n.y()*0.2, c.z() + n.z()*0.2);
             rl.drawLine3D(vc, ve, .blue);
        }

        // Draw current support point
        const sp = epa_s.current_support;
        rl.drawSphere(.init(sp.x(), sp.y(), sp.z()), 0.05, .red);
    }
}

fn drawUI(gjk_state: GJKStepper, epa_state: ?EPAState) void {
    rl.drawText("GJK/EPA Visualizer", 10, 10, 20, .black);
    if (gjk_state.finished) {
         if (gjk_state.intersect) {
             rl.drawText("GJK: INTERSECT", 10, 40, 20, .green);
         } else {
             rl.drawText("GJK: NO INTERSECT", 10, 40, 20, .red);
         }
    } else {
         rl.drawText("GJK: Running...", 10, 40, 20, .blue);
    }

    if (epa_state) |epa_s| {
         if (epa_s.finished) {
              rl.drawText("EPA: FINISHED", 10, 70, 20, .green);
         } else {
              var buf: [64]u8 = undefined;
              const s = std.fmt.bufPrintZ(&buf, "EPA: Step {} (V={})", .{epa_s.iter, epa_s.edges_count}) catch "";
              rl.drawText(s, 10, 70, 20, .blue);
         }
    }

    rl.drawText("Controls: P (Pause/Resume), N (Next Step), R (Reset)", 10, 420, 20, .black);
}
