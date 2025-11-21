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

    std.debug.print("Example complete.\n", .{});

    var app = try App.init("raylib [core] example - 3d camera free");
    defer app.deinit();

    try createScene(&app.world);

    while (!app.shouldClose()) {
        try app.update();

        app.beginDraw();
        app.drawScene();

        rl.drawRectangle(10, 10, 320, 93, .fade(.sky_blue, 0.5));
        rl.drawRectangleLines(10, 10, 320, 93, .blue);

        rl.drawText("Free camera default controls:", 20, 20, 10, .black);
        rl.drawText("- Mouse Wheel to Zoom in-out", 40, 40, 10, .dark_gray);
        rl.drawText("- Mouse Wheel Pressed to Pan", 40, 60, 10, .dark_gray);
        rl.drawText("- Z to zoom to (0, 0, 0)", 40, 80, 10, .dark_gray);

        app.endDraw();
    }
}

fn createScene(world: *zphys.World) !void {
    var ground = zphys.BodyDef.default();
    ground.shape = zphys.shape.newBox(math.vec3(5, 0.5, 5));
    ground.position = math.vec3(0, -0.5, 0);
    ground.inverseMass = 0.0;
    ground.friction = 0.9;
    ground.restitution = 0.2;
    _ = try world.createBody(ground);

    {
        var i: i32 = 0;
        while (i < 3) : (i += 1) {
            var d = zphys.BodyDef.default();
            d.shape = zphys.shape.newSphere(0.5);
            d.position = math.vec3(0, 3 + @as(f32, @floatFromInt(i)) * 1.1, 0);
            d.inverseMass = 1.0;
            d.friction = 0.4;
            d.restitution = 0.6;
            _ = try world.createBody(d);
        }
    }

    {
        var b1 = zphys.BodyDef.default();
        b1.shape = zphys.shape.newBox(math.vec3(0.5, 0.5, 0.5));
        b1.position = math.vec3(1.0, 4.0, 0.0);
        b1.inverseMass = 1.0;
        b1.friction = 0.6;
        b1.restitution = 0.0;
        _ = try world.createBody(b1);

        var b2 = zphys.BodyDef.default();
        b2.shape = zphys.shape.newBox(math.vec3(0.6, 0.4, 0.6));
        b2.position = math.vec3(1.8, 6.0, 0.1);
        b2.inverseMass = 1.0;
        b2.friction = 0.6;
        b2.restitution = 0.0;
        _ = try world.createBody(b2);
    }
}
