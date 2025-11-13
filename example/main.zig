const std = @import("std");
const math = @import("math");
const zphys = @import("zphys");
const rl = @import("raylib");
const DebugRenderer = @import("debug_renderer.zig").DebugRenderer;

pub fn main() !void {
    const a: math.Vec3 = math.vec3(1, 2, 3);
    const b: math.Vec3 = math.vec3(4, 5, 6);
    const c = a.add(&b);
    std.debug.print("math: a = {any}, b = {any}\n", .{ a, b });
    std.debug.print("math: a + b = {any}\n", .{ c });
    std.debug.print("math: dot(a, b) = {d}\n", .{ a.dot(&b) });

    const v: math.Vec2 = math.vec2(10, 20);
    std.debug.print("zphys: v = {any}, len = {d}\n", .{ v, v.len2() });

    std.debug.print("Example complete.\n", .{});

    const screenWidth = 800;
    const screenHeight = 450;

    rl.initWindow(screenWidth, screenHeight, "raylib [core] example - 3d camera free");
    defer rl.closeWindow();

    var camera = rl.Camera{
        .position = .init(10, 10, 10),
        .target = .init(0, 0, 0),
        .up = .init(0, 1, 0),
        .fovy = 45,
        .projection = .perspective,
    };

    var world = zphys.World.init(std.heap.page_allocator);
    defer world.deinit();

    var ground = zphys.BodyDef.default();
    ground.shape = zphys.shape.newBox(math.vec3(5, 0.5, 5));
    ground.position = math.vec3(0, -0.5, 0);
    ground.inverseMass = 0.0;
    ground.friction = 0.9;
    ground.restitution = 0.2;
    _ = try world.createBody(ground);

    //var i: i32 = 0;
    //while (i < 3) : (i += 1) {
    //    var d = zphys.BodyDef.default();
    //    d.shape = zphys.shape.newSphere(0.5);
    //    d.position = math.vec3(0, 3 + @as(f32, @floatFromInt(i)) * 1.1, 0);
    //    d.mass = 1.0;
    //    d.friction = 0.4;
    //    d.restitution = 0.6;
    //    _ = try world.createBody(d);
    //}

    {
        var b1 = zphys.BodyDef.default();
        b1.shape = zphys.shape.newBox(math.vec3(0.5, 0.5, 0.5));
        b1.position = math.vec3(1.0, 4.0, 0.0);
        b1.inverseMass = 1.0;
        b1.friction = 0.6;
        b1.restitution = 0.3;
        _ = try world.createBody(b1);

        //var b2 = zphys.BodyDef.default();
        //b2.shape = zphys.shape.newBox(math.vec3(0.6, 0.4, 0.6));
        //b2.position = math.vec3(1.8, 6.0, 0.1);
        //b2.inverseMass = 1.0;
        //b2.friction = 0.6;
        //b2.restitution = 0.3;
        //_ = try world.createBody(b2);
    }

    rl.disableCursor();
    rl.setTargetFPS(60);

    var paused: bool = false;
    var step_one: bool = false;

    const uvTex = try rl.loadTexture("example/resources/uvImageTexture.png");
    defer rl.unloadTexture(uvTex);

    var sky_tex: ?rl.Texture = null;
    defer if (sky_tex) |t| rl.unloadTexture(t);
    if (rl.loadImage("example/resources/sky.hdr")) |img| {
        if (rl.loadTextureFromImage(img)) |tex| {
            sky_tex = tex;
        } else |_| {}
        rl.unloadImage(img);
    } else |_| {}

    const albedo_index: usize = @intFromEnum(rl.MaterialMapIndex.albedo);
    const cube_mesh = rl.genMeshCube(1.0, 1.0, 1.0);
    var cube_model = try rl.loadModelFromMesh(cube_mesh);
    defer rl.unloadModel(cube_model);
    cube_model.materials[0].maps[albedo_index].texture = uvTex;

    const sphere_mesh = rl.genMeshSphere(1.0, 24, 24);
    var sphere_model = try rl.loadModelFromMesh(sphere_mesh);
    defer rl.unloadModel(sphere_model);
    sphere_model.materials[0].maps[albedo_index].texture = uvTex;

    while (!rl.windowShouldClose()) {
        camera.update(.free);

        if (rl.isKeyPressed(.space)) {
            paused = !paused;
        }
        
        if (rl.isKeyPressed(.right) and paused) {
            step_one = true;
        }

        if (!paused or step_one) {
            try world.step(1.0/60.0, 1);
            step_one = false;
        }

        if (rl.isKeyPressed(.z)) {
            camera.target = .init(0, 0, 0);
        }
        rl.beginDrawing();
        defer rl.endDrawing();

        rl.clearBackground(.ray_white);

        {
            if (sky_tex) |tex| {
                const src = rl.Rectangle.init(
                    0,
                    0,
                    @as(f32, @floatFromInt(tex.width)),
                    @as(f32, @floatFromInt(tex.height)),
                );
                const dst = rl.Rectangle.init(
                    0,
                    0,
                    @as(f32, @floatFromInt(screenWidth)),
                    @as(f32, @floatFromInt(screenHeight)),
                );
                rl.drawTexturePro(tex, src, dst, rl.Vector2.init(0, 0), 0.0, .white);
            } else {
                rl.drawRectangleGradientV(0, 0, screenWidth, screenHeight, .sky_blue, .ray_white);
            }
        }

        {
            camera.begin();
            defer camera.end();

            for (0..world.bodyCount()) |i| {
                const transform = world.getTransform(i);
                const shape = world.getShape(i);
                const trans_mat = math.Mat4x4.translate(transform.position);
                const rot_mat = math.Mat4x4.rotateByQuaternion(transform.orientation.normalize());
                switch (shape) {
                    .Box => |bx| {
                        const scale = bx.half_extents.mulScalar(2);
                        const scale_mat = math.Mat4x4.scale(scale);
                        const mat = trans_mat.mul(&rot_mat.mul(&scale_mat));
                        const rl_matrix = mathMat4ToRayLib(mat);
                        rl.drawMesh(cube_model.meshes[0], cube_model.materials[0], rl_matrix);
                    },
                    .Sphere => |sp| {
                        const scale_mat = math.Mat4x4.scale(math.vec3(sp.radius, sp.radius, sp.radius));
                        const mat = trans_mat.mul(&rot_mat.mul(&scale_mat));
                        const rl_matrix = mathMat4ToRayLib(mat);
                        rl.drawMesh(sphere_model.meshes[0], sphere_model.materials[0], rl_matrix);
                    },
                    .Line => |ln| {
                        const p1_local = ln.point_a.mulQuat(&transform.orientation);
                        const p2_local = ln.point_b.mulQuat(&transform.orientation);
                        const p1 = rl.Vector3.init(
                            transform.position.x() + p1_local.x(),
                            transform.position.y() + p1_local.y(),
                            transform.position.z() + p1_local.z(),
                        );
                        const p2 = rl.Vector3.init(
                            transform.position.x() + p2_local.x(),
                            transform.position.y() + p2_local.y(),
                            transform.position.z() + p2_local.z(),
                        );
                        rl.drawLine3D(p1, p2, .black);
                    },
                }
            }

            DebugRenderer.drawContacts(world.temp.contactSlice());
            DebugRenderer.drawManifolds(world.temp.manifoldSlice());

            rl.drawGrid(10, 1);
        }

        rl.drawRectangle(10, 10, 320, 93, .fade(.sky_blue, 0.5));
        rl.drawRectangleLines(10, 10, 320, 93, .blue);

        rl.drawText("Free camera default controls:", 20, 20, 10, .black);
        rl.drawText("- Mouse Wheel to Zoom in-out", 40, 40, 10, .dark_gray);
        rl.drawText("- Mouse Wheel Pressed to Pan", 40, 60, 10, .dark_gray);
        rl.drawText("- Z to zoom to (0, 0, 0)", 40, 80, 10, .dark_gray);
        
        DebugRenderer.drawDebugInfo(paused);
    }
}

fn mathMat4ToRayLib(math_mat: math.Mat4x4) rl.Matrix {
    return rl.Matrix{
        .m0 = math_mat.v[0].v[0],
        .m1 = math_mat.v[0].v[1],
        .m2 = math_mat.v[0].v[2],
        .m3 = math_mat.v[0].v[3],
        .m4 = math_mat.v[1].v[0],
        .m5 = math_mat.v[1].v[1],
        .m6 = math_mat.v[1].v[2],
        .m7 = math_mat.v[1].v[3],
        .m8 = math_mat.v[2].v[0],
        .m9 = math_mat.v[2].v[1],
        .m10 = math_mat.v[2].v[2],
        .m11 = math_mat.v[2].v[3],
        .m12 = math_mat.v[3].v[0],
        .m13 = math_mat.v[3].v[1],
        .m14 = math_mat.v[3].v[2],
        .m15 = math_mat.v[3].v[3],
    };
}
