pub const World = @import("world.zig").World;
pub const BodyDef = @import("body.zig").BodyDef;
pub const Contact = @import("collision/contact.zig");
pub const shape = @import("collision/shape.zig");
pub const collision = @import("collision/collision.zig");

const std = @import("std");

test {
    std.testing.refAllDeclsRecursive(@This());
}
