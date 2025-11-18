const std = @import("std");
const math = @import("math");
const MotionComp = @import("body.zig").MotionComp;
const contact = @import("collision/contact.zig");
const PenetrationConstraint = contact.PenetrationConstraint;


/// Solve Jacobian constraint(-n, -(r_0 x n, n, r_2 x n)
inline fn solveContactConstraint(constraint: *contact.PenetrationConstraint, motionA: *MotionComp, motionB: *MotionComp) void {
    var jv = constraint.n.dot(&motionB.velocity.sub(&motionA.velocity));
    jv -= constraint.r1.cross(&constraint.n).dot(&motionA.angularVelocity);
    jv += constraint.r2.cross(&constraint.n).dot(&motionB.angularVelocity);

    var impulse = -constraint.inverse_effective_mass * (jv - constraint.velocity_bias);
    {
        const new_accumulated = @max(0.0, constraint.accumulated_impulse + impulse);
        impulse = new_accumulated - constraint.accumulated_impulse;
        constraint.accumulated_impulse = new_accumulated;
    }

    // integrate velocity
    if (impulse != 0) {
        motionA.velocity = motionA.velocity.sub(&constraint.n.mulScalar(impulse * constraint.inv_mass_a));
        motionA.angularVelocity = motionA.angularVelocity.sub(&constraint.invert_inertia_n_x_r1.mulScalar(impulse));
        motionB.velocity = motionB.velocity.add(&constraint.n.mulScalar(impulse * constraint.inv_mass_b));
        motionB.angularVelocity = motionB.angularVelocity.add(&constraint.invert_inertia_n_x_r2.mulScalar(impulse));
    }

    // todo: Can we optimize the friction calculation by doing it in a single step?
    // todo: Remove code repetition
    const out_tangent1 = constraint.n.getNormalizePerpendicular();

    jv = out_tangent1.dot(&motionB.velocity.sub(&motionA.velocity));
    jv -= constraint.r1.cross(&out_tangent1).dot(&motionA.angularVelocity);
    jv += constraint.r2.cross(&out_tangent1).dot(&motionB.angularVelocity);

    impulse = -constraint.inverse_effective_mass * (jv - constraint.surface_velocity_bias1);
    {
        const new_accumulated = @max(0.0, constraint.accumulated_impulse + impulse);
        impulse = new_accumulated - constraint.accumulated_impulse;
        constraint.accumulated_impulse = new_accumulated;
    }

    // integrate velocity
    if (impulse != 0) {
        motionA.velocity = motionA.velocity.sub(&out_tangent1.mulScalar(impulse * constraint.inv_mass_a));
        motionA.angularVelocity = motionA.angularVelocity.sub(&constraint.invert_inertia_n_x_r1.mulScalar(impulse));
        motionB.velocity = motionB.velocity.add(&out_tangent1.mulScalar(impulse * constraint.inv_mass_b));
        motionB.angularVelocity = motionB.angularVelocity.add(&constraint.invert_inertia_n_x_r2.mulScalar(impulse));
    }

    const out_tangent2 = constraint.n.cross(&out_tangent1);
    jv = out_tangent2.dot(&motionB.velocity.sub(&motionA.velocity));
    jv -= constraint.r1.cross(&out_tangent2).dot(&motionA.angularVelocity);
    jv += constraint.r2.cross(&out_tangent2).dot(&motionB.angularVelocity);

    impulse = -constraint.inverse_effective_mass * (jv - constraint.surface_velocity_bias1);
    {
        const new_accumulated = @max(0.0, constraint.accumulated_impulse + impulse);
        impulse = new_accumulated - constraint.accumulated_impulse;
        constraint.accumulated_impulse = new_accumulated;
    }

    // integrate velocity
    if (impulse != 0) {
        motionA.velocity = motionA.velocity.sub(&out_tangent2.mulScalar(impulse * constraint.inv_mass_a));
        motionA.angularVelocity = motionA.angularVelocity.sub(&constraint.invert_inertia_n_x_r1.mulScalar(impulse));
        motionB.velocity = motionB.velocity.add(&out_tangent2.mulScalar(impulse * constraint.inv_mass_b));
        motionB.angularVelocity = motionB.angularVelocity.add(&constraint.invert_inertia_n_x_r2.mulScalar(impulse));
    }
}

/// Iteratively solve all penetration constraints
/// This is the main constraint solver that should be called from the physics step
pub fn solveConstraints(
    motion: []  MotionComp,
    constraints: []  PenetrationConstraint,
    iterations: u32
) void {
    var iteration: u32 = 0;
    while (iteration < iterations) : (iteration += 1) {
        for (constraints) |*constraint| {
            solveContactConstraint(
                constraint,
                &motion[constraint.body_a],
                &motion[constraint.body_b]
            );
        }
    }
}
