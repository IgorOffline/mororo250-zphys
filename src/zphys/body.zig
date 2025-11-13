const math = @import("math");
const shape = @import("collision/shape.zig");
const Shape = shape.Shape;

// Todo: Should I move this into it's own file?
pub const BodyDef = struct {
    angularVelocity: math.Vec3,
    centerOfMass: math.Vec3,
    orientation: math.Quat,

    velocity: math.Vec3,
    position: math.Vec3,

    inverseMass: f32,
    friction: f32,
    restitution: f32,

    shape: Shape,

    pub fn default() BodyDef {
        return BodyDef{
            .angularVelocity = math.vec3(0, 0, 0),
            .orientation = math.Quat.identity(),
            .velocity = math.vec3(0, 0, 0),
            .position = math.vec3(0, 0, 0),
            .inverseMass = 0.0,
            .centerOfMass = math.vec3(0, 0, 0),
            .friction = 0.5,
            .restitution = 0.5,
            // default placeholder shape (unit sphere)
            .shape = shape.newSphere(1.0),
        };
    }
};

pub const MotionComp = struct {
    angularVelocity: math.Vec3,
    velocity: math.Vec3,

    pub fn fromDef(def: BodyDef) MotionComp {
        return MotionComp{
            .angularVelocity = def.angularVelocity,
            .velocity = def.velocity,
        };
    }
};

pub const TransformComp = struct {
    orientation: math.Quat,
    position: math.Vec3,

    pub fn fromDef(def: BodyDef) TransformComp {
        return TransformComp{
            .orientation = def.orientation,
            .position = def.position,
        };
    }
};


pub const PhysicsPropsComp = struct {
    inverseInertia: math.Mat3x3, // Todo use vec3 as we only need the diagonal line(is this true for all objects?)
    centerOfMass: math.Vec3,
    inverseMass: f32,
    friction: f32,
    restitution: f32,

    pub fn fromDef(def: BodyDef, inverse_inertia: math.Mat3x3) PhysicsPropsComp {
        const inverse_mass = if (def.inverseMass == 0) 0 else 1.0 / def.inverseMass;
        return PhysicsPropsComp{
            .inverseInertia = inverse_inertia,
            .inverseMass = inverse_mass,
            .centerOfMass = def.centerOfMass,
            .friction = def.friction,
            .restitution = def.restitution,
        };
    }
};

pub const BodyComponents = struct {
    motion: MotionComp,
    transform: TransformComp,
    physics_props: PhysicsPropsComp,
    shape: Shape,
};

pub fn componentsFromDef(def: BodyDef) struct { MotionComp, TransformComp, PhysicsPropsComp, Shape } {
    const inverse_mass = if (def.inverseMass == 0) 0 else 1.0 / def.inverseMass;
    
    var inv_local = math.Mat3x3.init(&math.vec3(0, 0, 0), &math.vec3(0, 0, 0), &math.vec3(0, 0, 0));
    if (inverse_mass != 0) {
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
                const coef = def.inverseMass / 12.0;
                const xlength = 2 * b.half_extents.x();
                const ylength = 2 * b.half_extents.y();
                const zlength = 2 * b.half_extents.z();
                const Ixx: f32 = coef * (ylength * ylength + zlength * zlength);
                const Iyy: f32 = coef * (xlength * xlength + zlength * zlength);
                const Izz: f32 = coef * (xlength * xlength + ylength * ylength);
                const invIxx: f32 = if (Ixx > 0) 1.0 / Ixx else 0.0;
                const invIyy: f32 = if (Iyy > 0) 1.0 / Iyy else 0.0;
                const invIzz: f32 = if (Izz > 0) 1.0 / Izz else 0.0;
                inv_local = math.Mat3x3.init(
                    &math.vec3(invIxx, 0, 0),
                    &math.vec3(0, invIyy, 0),
                    &math.vec3(0, 0, invIzz),
                );
            },
            else => {},
        }
    }

    return .{
        MotionComp.fromDef(def),
        TransformComp.fromDef(def),
        PhysicsPropsComp.fromDef(def, inv_local),
        def.shape,
    };
}
