const std = @import("std");
const math = @import("math");
const zphys = @import("zphys");
const rl = @import("raylib");
const DebugRenderer = @import("debug_renderer.zig").DebugRenderer;
const SceneRenderer = @import("scene_renderer.zig").SceneRenderer;

pub fn main() !void {
    const screenWidth = 1280;
    const screenHeight = 720;

    rl.initWindow(screenWidth, screenHeight, "zphys - Pyramid Stack");
    defer rl.closeWindow();

    var camera = rl.Camera{
        .position = .init(10, 10, 20),
        .target = .init(0, 4, 0),
        .up = .init(0, 1, 0),
        .fovy = 45,
        .projection = .perspective,
    };

    var world = zphys.World.init(std.heap.page_allocator);
    defer world.deinit();

    // Ground
    var ground = zphys.BodyDef.default();
    ground.shape = zphys.shape.newBox(math.vec3(50, 1.0, 50));
    ground.position = math.vec3(0, -1.0, 0);
    ground.inverseMass = 0.0;
    ground.friction = 1.0;
    _ = try world.createBody(ground);

    // Pyramid
    const stack_height = 10;
    const box_size = 1.0;
    const gap = 0.05;

    var i: usize = 0;
    while (i < stack_height) : (i += 1) {
        const row_count = stack_height - i;
        const level_y = 0.5 + @as(f32, @floatFromInt(i)) * (box_size + gap);

        // Center the row along X
        const row_width = @as(f32, @floatFromInt(row_count)) * (box_size + gap) - gap;
        const start_x = -row_width * 0.5 + box_size * 0.5;

        var j: usize = 0;
        while (j < row_count) : (j += 1) {
            const x = start_x + @as(f32, @floatFromInt(j)) * (box_size + gap);
            
            var box = zphys.BodyDef.default();
            box.shape = zphys.shape.newBox(math.vec3(box_size * 0.5, box_size * 0.5, box_size * 0.5));
            box.position = math.vec3(x, level_y, 0);
            box.inverseMass = 1.0; // Dynamic
            box.friction = 0.5;
            box.restitution = 0.0;
            _ = try world.createBody(box);
        }
    }

    rl.disableCursor();
    rl.setTargetFPS(60);

    var paused: bool = false;
    var step_one: bool = false;

    var scene_renderer = try SceneRenderer.init();
    defer scene_renderer.deinit();

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
            camera.target = .init(0, 4, 0);
        }

        rl.beginDrawing();
        defer rl.endDrawing();

        SceneRenderer.drawSky();
        camera.begin();
        defer camera.end();

        scene_renderer.drawWorld(&world);

        // Draw contacts if needed, maybe too cluttered for pyramid
        // DebugRenderer.drawContacts(world.temp.contactSlice());
        
        rl.drawText("Pyramid Stack Test", 20, 20, 20, .black);
        DebugRenderer.drawDebugInfo(paused);
    }
}
