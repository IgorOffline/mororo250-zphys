const std = @import("std");
const builtin = @import("builtin");

pub fn main() !void {
    var mem: [4096]u8 = undefined;
    var fpa = std.heap.FixedBufferAllocator.init(&mem);
    const allocator = fpa.allocator();
    var it = try std.process.argsWithAllocator(allocator);
    defer it.deinit();

    var test_filter_opt: ?[]const u8 = null;
    while (it.next()) |arg| {
        if (std.mem.eql(u8, arg, "--test-filter")) {
            if (it.next()) |value| {
                test_filter_opt = try allocator.dupe(u8, value);
            } else {
                return error.MissingValueForTestFilter;
            }
        } else if (std.mem.startsWith(u8, arg, "--test-filter=")) {
            test_filter_opt = arg["--test-filter=".len..];
        } else {
            // handle other args
        }
    }

    var stderr_buf: [1024]u8 = undefined;
    var stderr_writter = std.fs.File.stderr().writer(&stderr_buf);
    var stderr = &stderr_writter.interface;

    const tty_config = std.Io.tty.detectConfig(.stderr());

    var total_count: usize = 0;
    for (builtin.test_functions) |t| {
        if (test_filter_opt) |test_filter| {
            if (std.mem.indexOf(u8, t.name, test_filter) == null)
                continue;
        }
        total_count += 1;
    }

    var pass_count: usize = 0;
    var fail_count: usize = 0;
    var leak_count: usize = 0;
    var total_elapsed: i64 = 0;

    for (builtin.test_functions) |t| {
        if (test_filter_opt) |test_filter| {
            if (std.mem.indexOf(u8, t.name, test_filter) == null)
                continue;
        }

        std.testing.allocator_instance = .{};
        const name = extractName(t);

        // Capture stderr output during test ex

        const start = std.time.milliTimestamp();
        const result = t.func();
        const elapsed = std.time.milliTimestamp() - start;
        total_elapsed += elapsed;
        if (std.testing.allocator_instance.deinit() == .leak) {
            try printRed(&tty_config, stderr, "{s} leaked memory\n", .{name});
            leak_count += 1;
        }
        if (result) |_| {
            pass_count += 1;
            try printGreen(&tty_config, stderr, "{s} - ({d}ms)\n", .{ name, elapsed });
        } else |err| {
            fail_count += 1;
            try printRed(&tty_config, stderr, "{s} - {}\n", .{ name, err });
        }
    }
    try stderr.print("\nSummary: total={d}, passed={d}, failed={d}, leaked={d}, time={d}ms\n", .{ total_count, pass_count, fail_count, leak_count, total_elapsed });
    try stderr.flush();
}

fn extractName(t: std.builtin.TestFn) []const u8 {
    const marker = std.mem.lastIndexOf(u8, t.name, ".test.") orelse return t.name;
    return t.name[marker + 6 ..];
}

fn printGreen(config: *const std.Io.tty.Config, writer: *std.Io.Writer, comptime message: []const u8, args: anytype) !void {
    try config.setColor(writer, .green);
    try writer.print(message, args);
    try config.setColor(writer, .white);
    try writer.flush();
}

fn printRed(config: *const std.Io.tty.Config, writer: *std.Io.Writer, comptime message: []const u8, args: anytype) !void {
    try config.setColor(writer, .red);
    try writer.print(message, args);
    try config.setColor(writer, .white);
    try writer.flush();
}
