const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const exe = b.addExecutable(.{
        .name = "main",
        .target = target,
        .optimize = optimize,
    });

    exe.addIncludePath(b.path("external/nuklear"));
    exe.addIncludePath(b.path("external/stb"));
    exe.addCSourceFile(.{ .file = b.path("main.c"), .flags = &.{} });

    const sdl_dep = b.dependency("sdl", .{
        .target = target,
        .optimize = optimize,
        //.preferred_link_mode = .static, // or .dynamic
    });
    const sdl_lib = sdl_dep.artifact("SDL3");
    exe.linkLibrary(sdl_lib);

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    b.step("run", "Run the application").dependOn(&run_cmd.step);
}
