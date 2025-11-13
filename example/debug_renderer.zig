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
    
    pub fn drawContacts(contacts: []const zphys.Contact.Contact) void {
        for (contacts) |contact| {
            const point_a_rl = rl.Vector3.init(contact.point_a.x(), contact.point_a.y(), contact.point_a.z());
            const point_b_rl = rl.Vector3.init(contact.point_b.x(), contact.point_b.y(), contact.point_b.z());
            
            rl.drawSphere(point_a_rl, 0.05, .blue);
            rl.drawSphere(point_b_rl, 0.05, .red);
            
            const normal_end = contact.point_a.add(&contact.normal.mulScalar(0.5));
            const normal_end_rl = rl.Vector3.init(normal_end.x(), normal_end.y(), normal_end.z());
            rl.drawLine3D(point_a_rl, normal_end_rl, .yellow);
        }
    }
    
    pub fn drawManifolds(manifolds: []const zphys.Contact.ContactManifold) void {
        for (manifolds) |manifold| {
            var world_points_a: [4]rl.Vector3 = undefined;
            var world_points_b: [4]rl.Vector3 = undefined;

            var i: u32 = 0;
            while (i < manifold.length) : (i += 1) {
                world_points_a[i] = rl.Vector3.init(manifold.contact_points_a[i].x(), manifold.contact_points_a[i].y(), manifold.contact_points_a[i].z());
                world_points_b[i] = rl.Vector3.init(manifold.contact_points_b[i].x(), manifold.contact_points_b[i].y(), manifold.contact_points_b[i].z());

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
