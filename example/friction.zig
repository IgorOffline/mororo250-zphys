const std = @import("std");
const math = @import("math");
const zphys = @import("zphys");
const rl = @import("raylib");
const App = @import("shared/app.zig").App;

pub fn main() !void {
    const a: math.Vec3 = math.vec3(1, 2, 3);
    const b: math.Vec3 = math.vec3(4, 5, 6);
    const c = a.add(&b);
    std.debug.print("math: a = {any}, b = {any}\n", .{ a, b });
    std.debug.print("math: a + b = {any}\n", .{ c });
    std.debug.print("math: dot(a, b) = {d}\n", .{ a.dot(&b) });

    const v: math.Vec2 = math.vec2(10, 20);
    std.debug.print("zphys: v = {any}, len = {d}\n", .{ v, v.len2() });

    std.debug.print("Friction test example.\n", .{});

    var app = try App.init("zphys - Friction Test");
    defer app.deinit();

    try createScene(&app.world);

    while (!app.shouldClose()) {
        try app.update();

        app.beginDraw();
        app.drawScene();

        rl.drawRectangle(10, 10, 320, 133, .fade(.sky_blue, 0.5));
        rl.drawRectangleLines(10, 10, 320, 133, .blue);

        rl.drawText("Friction Test Scene", 20, 20, 10, .black);
        rl.drawText("- Spheres (friction 0.0-0.8)", 40, 40, 10, .dark_gray);
        rl.drawText("- Boxes (friction 0.0-0.8)", 40, 60, 10, .dark_gray);
        rl.drawText("- Mouse Wheel to Zoom", 40, 80, 10, .dark_gray);
        rl.drawText("- Mouse Wheel Pressed to Pan", 40, 100, 10, .dark_gray);
        rl.drawText("- Z to zoom to (0, 0, 0)", 40, 120, 10, .dark_gray);

        app.endDraw();
    }
}

fn createScene(world: *zphys.World) !void {
    try createRamp(world);
    try createFrictionObjects(world);
    try createContainer(world);
}

fn createRamp(world: *zphys.World) !void {
    // Create ramp - a long box rotated at an angle
    const ramp_angle = std.math.pi / 6.0; // 30 degrees
    var ramp = zphys.BodyDef.default();
    ramp.shape = zphys.shape.newBox(math.vec3(10.0, 0.5, 4.0));
    ramp.position = math.vec3(0, 0, 0);
    ramp.orientation = math.Quat.fromAxisAngle(math.vec3(0, 0, 1), ramp_angle);
    ramp.inverseMass = 0.0;
    ramp.friction = 0.5;
    ramp.restitution = 0.0;
    _ = try world.createBody(ramp);
}

fn createFrictionObjects(world: *zphys.World) !void {
    // Friction test objects - Spheres and boxes with different friction values
    const friction_values = [_]f32{ 0.0, 0.2, 0.4, 0.6, 0.8 };

    // Calculate the top of the ramp position
    const ramp_top_x = 7.0;
    const ramp_top_y = 8.0;
    const ramp_angle = std.math.pi / 6.0;

    // Create spheres and boxes side by side at the top of the ramp
    for (friction_values, 0..) |friction, i| {
        const z_offset = @as(f32, @floatFromInt(i)) * 1.5 - 3.0;

        // Create sphere
        var sphere = zphys.BodyDef.default();
        sphere.shape = zphys.shape.newSphere(0.4);
        sphere.position = math.vec3(ramp_top_x, ramp_top_y, z_offset - 0.5);
        sphere.inverseMass = 1.0;
        sphere.friction = friction;
        sphere.restitution = 0.0;
        _ = try world.createBody(sphere);

        // Create box
        var box = zphys.BodyDef.default();
        box.shape = zphys.shape.newBox(math.vec3(0.4, 0.4, 0.4));
        box.position = math.vec3(ramp_top_x, ramp_top_y, z_offset + 0.5);
        box.orientation = math.Quat.fromAxisAngle(math.vec3(0, 0, 1), ramp_angle);
        box.inverseMass = 1.0;
        box.friction = friction;
        box.restitution = 0.1;
        _ = try world.createBody(box);
    }
}

fn createContainer(world: *zphys.World) !void {
    // Create container at the bottom to catch objects
    const container_pos = math.vec3(-12.0, -10.0, 0.0);
    const wall_thickness = 0.5;
    const half_width = 6.0;
    const half_length = 5.0;
    const wall_height = 3.0;

    // Floor
    var floor = zphys.BodyDef.default();
    floor.shape = zphys.shape.newBox(math.vec3(half_length, wall_thickness, half_width));
    floor.position = container_pos;
    floor.inverseMass = 0.0;
    floor.friction = 0.5;
    _ = try world.createBody(floor);

    // Back Wall (Stopping wall)
    var back_wall = zphys.BodyDef.default();
    back_wall.shape = zphys.shape.newBox(math.vec3(wall_thickness, wall_height, half_width));
    back_wall.position = container_pos.add(&math.vec3(-half_length - wall_thickness, wall_height, 0));
    back_wall.inverseMass = 0.0;
    _ = try world.createBody(back_wall);

    // Front Wall (Near ramp) - Make it lower so they fall in more easily?
    // Or assume ramp ends above it.
    // Ramp end: (-8.6, -5). Container floor: -10.
    // So they drop 5 units.
    var front_wall = zphys.BodyDef.default();
    front_wall.shape = zphys.shape.newBox(math.vec3(wall_thickness, wall_height, half_width));
    front_wall.position = container_pos.add(&math.vec3(half_length + wall_thickness, wall_height, 0));
    front_wall.inverseMass = 0.0;
    _ = try world.createBody(front_wall);

    // Side Wall 1
    var side_wall1 = zphys.BodyDef.default();
    side_wall1.shape = zphys.shape.newBox(math.vec3(half_length, wall_height, wall_thickness));
    side_wall1.position = container_pos.add(&math.vec3(0, wall_height, half_width + wall_thickness));
    side_wall1.inverseMass = 0.0;
    _ = try world.createBody(side_wall1);

    // Side Wall 2
    var side_wall2 = zphys.BodyDef.default();
    side_wall2.shape = zphys.shape.newBox(math.vec3(half_length, wall_height, wall_thickness));
    side_wall2.position = container_pos.add(&math.vec3(0, wall_height, -half_width - wall_thickness));
    side_wall2.inverseMass = 0.0;
    _ = try world.createBody(side_wall2);
}
