const std = @import("std");
const rl = @import("raylib");
const zphys = @import("zphys");
const math = @import("math");
const raylibUtils = @import("raylibUtils.zig");

pub const SceneRenderer = struct {
    shader: rl.Shader,
    colSecondaryLoc: i32,
    cubeModel: rl.Model,
    sphereModel: rl.Model,

    const zigYellow = rl.Color.init(247, 164, 29, 255);
    const skyColor = rl.Color.sky_blue;
    
    const staticPrimary = rl.Color.init(200, 200, 200, 255);
    const staticSecondary = rl.Color.init(100, 100, 100, 255);

    pub fn init() !SceneRenderer {
        var self: SceneRenderer = undefined;

        // Load Shader
        // Note: raylib loadShader takes optional vsFileName and fsFileName.
        // We use default vertex shader (null) and our custom fragment shader.
        self.shader = try rl.loadShader(null, "example/resources/checker_shader.fs");
        self.colSecondaryLoc = rl.getShaderLocation(self.shader, "colSecondary");

        // Create Models
        const albedo_index: usize = @intFromEnum(rl.MaterialMapIndex.albedo);

        // Cube
        const cube_mesh = rl.genMeshCube(1.0, 1.0, 1.0);
        self.cubeModel = try rl.loadModelFromMesh(cube_mesh);
        self.cubeModel.materials[0].shader = self.shader;
        // Ensure white base color so shader modulation works
        self.cubeModel.materials[0].maps[albedo_index].color = rl.Color.white;

        // Sphere
        const sphere_mesh = rl.genMeshSphere(1.0, 24, 24);
        self.sphereModel = try rl.loadModelFromMesh(sphere_mesh);
        self.sphereModel.materials[0].shader = self.shader;
        self.sphereModel.materials[0].maps[albedo_index].color = rl.Color.white;

        return self;
    }

    pub fn deinit(self: *SceneRenderer) void {
        rl.unloadShader(self.shader);
        rl.unloadModel(self.cubeModel);
        rl.unloadModel(self.sphereModel);
    }

    pub fn drawSky() void {
        rl.clearBackground(skyColor);
    }

    pub fn drawWorld(self: SceneRenderer, world: *zphys.World) void {
        for (0..world.bodyCount()) |i| {
            const transform = world.getTransform(i);
            const shape = world.getShape(i);
            const props = world.getPhysicsProps(i);
            const trans_mat = math.Mat4x4.translate(transform.position);
            const rot_mat = math.Mat4x4.rotateByQuaternion(transform.orientation.normalize());

            var colorPrimary: rl.Color = undefined;
            var colorSecondary: rl.Color = undefined;
            
            if (props.inverseMass == 0) {
                colorPrimary = staticPrimary;
                colorSecondary = staticSecondary;
            } else {
                colorPrimary = rl.Color.white;
                colorSecondary = zigYellow;
            }
            
            const colSecFloats = raylibUtils.colorToFloatArray(colorSecondary);
            rl.setShaderValue(self.shader, self.colSecondaryLoc, &colSecFloats, .vec4);
            const albedo_index: usize = @intFromEnum(rl.MaterialMapIndex.albedo);

            switch (shape) {
                .Box => |bx| {
                    const scale = bx.half_extents.mulScalar(2);
                    const scale_mat = math.Mat4x4.scale(scale);
                    const mat = trans_mat.mul(&rot_mat.mul(&scale_mat));
                    const rl_matrix = raylibUtils.mathMat4ToRayLib(mat);
                    
                    self.cubeModel.materials[0].maps[albedo_index].color = colorPrimary;
                    rl.drawMesh(self.cubeModel.meshes[0], self.cubeModel.materials[0], rl_matrix);
                },
                .Sphere => |sp| {
                    const scale_mat = math.Mat4x4.scale(math.vec3(sp.radius, sp.radius, sp.radius));
                    const mat = trans_mat.mul(&rot_mat.mul(&scale_mat));
                    const rl_matrix = raylibUtils.mathMat4ToRayLib(mat);

                    self.sphereModel.materials[0].maps[albedo_index].color = colorPrimary;
                    rl.drawMesh(self.sphereModel.meshes[0], self.sphereModel.materials[0], rl_matrix);
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
    }
};
