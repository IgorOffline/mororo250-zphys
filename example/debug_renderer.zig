const std = @import("std");
const math = @import("math");
const zphys = @import("zphys");
const rl = @import("raylib");

pub const DebugRenderer = struct {

    fn drawThickLine(start: rl.Vector3, end: rl.Vector3, thickness: f32, color: rl.Color) void {
        const direction = rl.Vector3.subtract(end, start);
        const length = rl.Vector3.length(direction);
        if (length < 0.001) return;
        
        rl.drawCylinderEx(start, end, thickness, thickness, 6, color);
    }
    
    pub fn drawContacts(contacts: []const zphys.Contact.Contact, bodies: []const zphys.Body) void {
        for (contacts) |contact| {
            if (contact.body_a >= bodies.len or contact.body_b >= bodies.len) continue;
            
            const body_a = bodies[contact.body_a];
            const body_b = bodies[contact.body_b];
            
            // Transform local points to world space
            const world_point_a_local = contact.point_local_a.mulQuat(&body_a.orientation);
            const world_point_a = body_a.position.add(&world_point_a_local);
            
            const world_point_b_local = contact.point_local_b.mulQuat(&body_b.orientation);
            const world_point_b = body_b.position.add(&world_point_b_local);
            
            // Draw contact points
            const point_a_rl = rl.Vector3.init(world_point_a.x(), world_point_a.y(), world_point_a.z());
            const point_b_rl = rl.Vector3.init(world_point_b.x(), world_point_b.y(), world_point_b.z());
            
            // Draw spheres at contact points - red for body A, green for body B
            rl.drawSphere(point_a_rl, 0.05, .red);
            rl.drawSphere(point_b_rl, 0.05, .green);
            
            // Draw contact normal
            const normal_end = world_point_a.add(&contact.normal.mulScalar(0.5));
            const normal_end_rl = rl.Vector3.init(normal_end.x(), normal_end.y(), normal_end.z());
            rl.drawLine3D(point_a_rl, normal_end_rl, .yellow);
        }
    }
    
    pub fn drawManifolds(manifolds: []const zphys.Contact.ContactManifold, bodies: []const zphys.Body) void {
        for (manifolds) |manifold| {
            if (manifold.body_a >= bodies.len or manifold.body_b >= bodies.len) continue;
            
            const body_a = bodies[manifold.body_a];
            const body_b = bodies[manifold.body_b];
            
            var world_points_a: [4]rl.Vector3 = undefined;
            var world_points_b: [4]rl.Vector3 = undefined;
            
            std.debug.print("\n=== RENDERING MANIFOLD ===\n", .{});
            std.debug.print("Body A ID: {}, Body B ID: {}\n", .{manifold.body_a, manifold.body_b});
            std.debug.print("Body A pos: ({d:.3}, {d:.3}, {d:.3})\n", .{body_a.position.x(), body_a.position.y(), body_a.position.z()});
            std.debug.print("Body B pos: ({d:.3}, {d:.3}, {d:.3})\n", .{body_b.position.x(), body_b.position.y(), body_b.position.z()});
            
            var i: u32 = 0;
            while (i < manifold.length) : (i += 1) {
                std.debug.print("\n  Converting point {}:\n", .{i});
                std.debug.print("    Stored Local A: ({d:.3}, {d:.3}, {d:.3})\n", .{
                    manifold.contact_points_a[i].x(),
                    manifold.contact_points_a[i].y(),
                    manifold.contact_points_a[i].z()
                });
                
                const point_a_local = manifold.contact_points_a[i].mulQuat(&body_a.orientation);
                std.debug.print("    After Quat rotation: ({d:.3}, {d:.3}, {d:.3})\n", .{
                    point_a_local.x(),
                    point_a_local.y(),
                    point_a_local.z()
                });
                
                const world_point_a = body_a.position.add(&point_a_local);
                std.debug.print("    Final World A: ({d:.3}, {d:.3}, {d:.3})\n", .{
                    world_point_a.x(),
                    world_point_a.y(),
                    world_point_a.z()
                });
                
                const point_b_local = manifold.contact_points_b[i].mulQuat(&body_b.orientation);
                const world_point_b = body_b.position.add(&point_b_local);
                std.debug.print("    Final World B: ({d:.3}, {d:.3}, {d:.3})\n", .{
                    world_point_b.x(),
                    world_point_b.y(),
                    world_point_b.z()
                });
                
                world_points_a[i] = rl.Vector3.init(world_point_a.x(), world_point_a.y(), world_point_a.z());
                world_points_b[i] = rl.Vector3.init(world_point_b.x(), world_point_b.y(), world_point_b.z());
                
                rl.drawSphere(world_points_a[i], 0.08, .blue);
                rl.drawSphere(world_points_b[i], 0.08, .red);
                
                drawThickLine(world_points_a[i], world_points_b[i], 0.02, .purple);
            }
            
            if (manifold.length >= 2) {
                var j: u32 = 0;
                while (j < manifold.length) : (j += 1) {
                    const next = (j + 1) % manifold.length;
                    drawThickLine(world_points_a[j], world_points_a[next], 0.025, .green);
                    drawThickLine(world_points_b[j], world_points_b[next], 0.025, .orange);
                }
            }
        }
    }
    
    pub fn drawDebugInfo(paused: bool) void {
        const y_offset: i32 = 120;
        rl.drawRectangle(10, y_offset, 320, 73, .fade(.lime, 0.5));
        rl.drawRectangleLines(10, y_offset, 320, 73, .dark_green);
        
        rl.drawText("Debug Controls:", 20, y_offset + 10, 10, .black);
        rl.drawText("- SPACE: Toggle Pause", 40, y_offset + 25, 10, .dark_gray);
        rl.drawText("- RIGHT: Step One Frame (when paused)", 40, y_offset + 40, 10, .dark_gray);
        
        const status = if (paused) "PAUSED" else "RUNNING";
        const status_color = if (paused) rl.Color.red else rl.Color.green;
        rl.drawText(status, 40, y_offset + 55, 10, status_color);
    }
};
