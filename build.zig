const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const math_mod = b.addModule("math", .{
        .root_source_file = b.path("src/math/math.zig"),
        .target = target,
        .optimize = optimize,
    });

    const physics_mod = b.addModule("zphys", .{
        .root_source_file = b.path("src/zphys/zphys.zig"),
        .target = target,
        .optimize = optimize,
    });
    // allow the zphys module to import the math module via @import("math")
    physics_mod.addImport("math", math_mod);

    const raylib_dep = b.dependency("raylib_zig", .{
        .target = target,
        .optimize = optimize,
    });
    const raylib = raylib_dep.module("raylib");
    const raylib_artifact = raylib_dep.artifact("raylib");

    var examples = [_]Example{
        .{ .name = "zphys_basic", .source = "example/basic.zig", .option = "basic" },
        .{ .name = "zphys_friction", .source = "example/friction.zig", .option = "friction" },
        .{ .name = "zphys_pyramid", .source = "example/pyramid.zig", .option = "pyramid" },
        .{ .name = "zphys_visualizer", .source = "example/gjk_epa_visualizer.zig", .option = "visualizer" },
    };

    // Create all example executables
    for (&examples) |*example| {
        example.exe = addExample(b, example.name, example.source, math_mod, physics_mod, raylib, raylib_artifact, target, optimize);
    }

    const test_step = b.step("test", "Run tests");

    // tests for the math library
    const math_mod_tests = b.addTest(.{
        .name = "math_mod_tests",
        .root_module = math_mod,
        .test_runner = .{.path = b.path("src/test_runner.zig"), .mode = .simple },
    });
    const run_math_tests = b.addRunArtifact(math_mod_tests);
    test_step.dependOn(&run_math_tests.step);

    // tests for the physics library
    const physics_mod_tests = b.addTest(.{
        .name = "physics_mod_tests",
        .root_module = physics_mod,
        .test_runner = .{.path = b.path("src/test_runner.zig"), .mode = .simple },
    });
    const run_physics_tests = b.addRunArtifact(physics_mod_tests);
    test_step.dependOn(&run_physics_tests.step);

    // Build options to select which example to run
    const run_step = b.step("run", "Run examples (use -Dbasic or -Dfriction)");
    
    var selected_example: ?*std.Build.Step.Compile = null;
    
    // Check which example option was selected
    for (&examples) |*example| {
        const is_selected = b.option(bool, example.option, b.fmt("Run {s} example", .{example.option})) orelse false;
        if (is_selected) {
            selected_example = example.exe;
            break;
        }
    }
    
    // If no option selected, default to basic example
    if (selected_example == null) {
        selected_example = examples[0].exe;
    }
    
    // Run the selected example
    if (selected_example) |exe| {
        const run_cmd = b.addRunArtifact(exe);
        run_cmd.step.dependOn(b.getInstallStep());
        if (b.args) |args| {
            run_cmd.addArgs(args);
        }
        run_step.dependOn(&run_cmd.step);
    }

    if (b.args) |args| {
        run_math_tests.addArgs(args);
        run_physics_tests.addArgs(args);
    }

    const check = examples[0].exe;

    // This is a test of making zls work in tests.
    const check_test = b.addTest(.{
        .root_module = math_mod,
    });

    const check_step = b.step("check", "Check for zls analysis");
    check_step.dependOn(&check.step);
    check_step.dependOn(&check_test.step);
}

// Define examples structure
const Example = struct {
    name: []const u8,
    source: []const u8,
    option: []const u8,
    exe: *std.Build.Step.Compile = undefined,
};

// Helper function to add examples
fn addExample(
    builder: *std.Build,
    name: []const u8,
    source_file: []const u8,
    math: *std.Build.Module,
    physics: *std.Build.Module,
    rl: *std.Build.Module,
    rl_artifact: *std.Build.Step.Compile,
    tgt: std.Build.ResolvedTarget,
    opt: std.builtin.OptimizeMode,
) *std.Build.Step.Compile {
    const mod = builder.createModule(.{
        .root_source_file = builder.path(source_file),
        .target = tgt,
        .optimize = opt,
        .imports = &.{
            .{ .name = "math", .module = math },
            .{ .name = "zphys", .module = physics },
            .{ .name = "raylib", .module = rl },
        },
    });

    const exe = builder.addExecutable(.{
        .name = name,
        .root_module = mod,
    });
    exe.linkLibrary(rl_artifact);
    builder.installArtifact(exe);
    return exe;
}
