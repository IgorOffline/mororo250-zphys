const std = @import("std");
const math = @import("math");
const Body = @import("body.zig").Body;
const BodyDef = @import("body.zig").BodyDef;
const collision = @import("collision/collision.zig");

// Todo: Separate static bodies into different list for optimization reasons
pub const World = struct {
    allocator: std.mem.Allocator,
    bodies: std.ArrayList(Body),
    gravity: math.Vec3,
    temp: WorldTemp,

    pub fn init(allocator: std.mem.Allocator) World {
        return World.initWithGravity(allocator, math.vec3(0, -9.81, 0));
    }

    pub fn initWithGravity(allocator: std.mem.Allocator, gravity: math.Vec3) World {
        return .{
            .allocator = allocator,
            .bodies = .{},
            .gravity = gravity,
            .temp = WorldTemp.init(allocator),
        };
    }

    pub fn deinit(self: *World) void {
        self.temp.deinit();
        self.bodies.deinit(self.allocator);
    }

    pub fn createBody(self: *World, def: BodyDef) !u32 {
        const id: u32 = @intCast(self.bodies.items.len);

        var body = Body.fromDef(def);

        var inv_local = math.Mat3x3.init(&math.vec3(0, 0, 0), &math.vec3(0, 0, 0), &math.vec3(0, 0, 0));
        if (body.inverseMass != 0) {
            switch (def.shape) {
                .Sphere => |s| {
                    const r2: f32 = s.radius * s.radius;
                    const I: f32 = (2.0 / 5.0) * def.inverseMass * r2;
                    const invI: f32 = if (I > 0) 1.0 / I else 0.0;
                    inv_local = math.Mat3x3.init(
                        &math.vec3(invI, 0, 0),
                        &math.vec3(0, invI, 0),
                        &math.vec3(0, 0, invI),
                    );
                },
                .Box => |b| {
                    const hx = b.half_extents.x();
                    const hy = b.half_extents.y();
                    const hz = b.half_extents.z();
                    const Ixx: f32 = (1.0 / 3.0) * def.inverseMass * (hy * hy + hz * hz);
                    const Iyy: f32 = (1.0 / 3.0) * def.inverseMass * (hx * hx + hz * hz);
                    const Izz: f32 = (1.0 / 3.0) * def.inverseMass * (hx * hx + hy * hy);
                    const invIxx: f32 = if (Ixx > 0) 1.0 / Ixx else 0.0;
                    const invIyy: f32 = if (Iyy > 0) 1.0 / Iyy else 0.0;
                    const invIzz: f32 = if (Izz > 0) 1.0 / Izz else 0.0;
                    inv_local = math.Mat3x3.init(
                        &math.vec3(invIxx, 0, 0),
                        &math.vec3(0, invIyy, 0),
                        &math.vec3(0, 0, invIzz),
                    );
                },
                else => {
                    unreachable;
                },
            }
        }

        body.inverseInertia = inv_local;


        try self.bodies.append(self.allocator, body);
        return id;
    }

    pub fn step(self: *World, timestep: f32, substep: u16) !void {
        std.debug.assert(substep > 0);
        const dt: f32 = timestep / @as(f32, @floatFromInt(substep));

        var substep_index: u16 = 0;
        try self.temp.ensureCapacity(self.bodies.items.len);
        while (substep_index < substep) : (substep_index += 1) {
            applyGravity(self, dt);

            self.temp.clear();
            collision.generateContacts(self.bodies.items, &self.temp.contacts);
            collision.solveVelocity(self.bodies.items, self.temp.contactSlice(), 10);

            integratePositions(self, dt);

            var iteration: u8 = 0;
            while (iteration < 10) : (iteration += 1) {
                collision.solvePosition(self.bodies.items, self.temp.contactSlice());
            }
        }
    }

    fn applyGravity(self: *World, dt: f32) void {
        var body_index: u16 = 0;
        while (body_index < self.bodies.items.len) : (body_index += 1) {
            var body = &self.bodies.items[body_index];
            if (body.inverseMass == 0) continue; // static
                const gravity_delta_velocity = self.gravity.mulScalar(dt);
            body.velocity = body.velocity.add(&gravity_delta_velocity);
        }
    }

    fn integratePositions(self: *World, dt: f32) void {
        for (0..self.bodies.items.len) |body_index| {
            var body = &self.bodies.items[body_index];
            if (body.inverseMass == 0) continue;
            const position_delta = body.velocity.mulScalar(dt);
            body.position = body.position.add(&position_delta);

            const omega = body.angularVelocity;
            const omega_len2 = omega.len2();
            if (omega_len2 > 1e-12) {
                const omega_len = std.math.sqrt(omega_len2);
                const axis = omega.mulScalar(1.0 / omega_len);
                const angle = omega_len * dt;
                const dq = math.Quat.fromAxisAngle(axis, angle);
                body.orientation = math.Quat.mul(&dq, &body.orientation).normalize();
            }
        }
    }
};

pub const WorldTemp = struct {
    allocator: std.mem.Allocator,
    contacts: std.ArrayList(collision.Contact),

    pub fn init(allocator: std.mem.Allocator) WorldTemp {
        return .{
            .allocator = allocator,
            .contacts = .{},
        };
    }

    pub fn deinit(self: *WorldTemp) void {
        self.contacts.deinit(self.allocator);
    }

    pub fn clear(self: *WorldTemp) void {
        self.contacts.clearRetainingCapacity();
    }

    pub fn ensureCapacity(self: *WorldTemp, bodies_count: usize) !void {
        if (bodies_count <= 1) return;
        const max_pairs = bodies_count * (bodies_count - 1) / 2;
        try self.contacts.ensureTotalCapacity(self.allocator, max_pairs);
    }

    pub fn contactSlice(self: *WorldTemp) []const collision.Contact {
        return self.contacts.items;
    }
};
