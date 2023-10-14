const std = @import("std");

const os_tag = @import("builtin").os.tag;

pub fn build(b: *std.Build) !void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const sdk_path = b.option([]const u8, "playdate-sdk-path", "Path to installed Playdate SDK");

    const playdate_mod = b.createModule(.{ .source_file = .{ .path = "src/playdate.zig" } });
    const modem_mod = b.addModule("modem", .{ .source_file = .{ .path = "src/modem.zig" } });

    const portaudio_dep = b.dependency(
        "portaudio",
        .{ .target = target, .optimize = optimize, .shared = false },
    );
    const portaudio_lib = portaudio_dep.artifact("portaudio");

    const example_send_game = try addPlaydateGameExe(b, .{ .path = "examples/playdate-send/main.zig" }, .{
        .optimize = optimize,
        .assets_dir = "examples/playdate-send/assets",
        .modules = &.{
            .{ .name = "modem", .module = modem_mod },
            .{ .name = "playdate", .module = playdate_mod },
        },
    });
    try addPlaydateSimulatorRun(b, "example-playdate-send", example_send_game, sdk_path);

    const log_level = b.option(std.log.Level, "log-level", "Log verbosity threshold") orelse .info;
    const options = b.addOptions();
    options.addOption(std.log.Level, "log_level", log_level);

    const pc_examples = .{
        .{ .name = "demodulator-dump", .src = "demodulator_dump.zig", .pa = true },
        .{ .name = "modulator-write-raw", .src = "modulator_write_raw.zig", .pa = false },
        .{ .name = "visualize-goertzel", .src = "visualize_goertzel.zig", .pa = true },
    };
    inline for (pc_examples) |ex| {
        const exe = b.addExecutable(.{
            .name = ex.name,
            .root_source_file = .{ .path = "examples/" ++ ex.src },
            .target = target,
            .optimize = optimize,
        });
        exe.addModule("modem", modem_mod);
        exe.addOptions("options", options);
        if (ex.pa) exe.linkLibrary(portaudio_lib);
        b.installArtifact(exe);

        const run = b.addRunArtifact(exe);
        run.step.dependOn(b.getInstallStep());
        if (b.args) |args| {
            run.addArgs(args);
        }

        const run_step = b.step("example-" ++ ex.name, "Run example " ++ ex.name);
        run_step.dependOn(&run.step);
    }

    const tests = b.addTest(.{
        .root_source_file = .{ .path = "src/modem.zig" },
        .target = target,
        .optimize = optimize,
    });
    tests.filter = b.option([]const u8, "test-filter", "Filter tests by name");

    const tests_run_cmd = b.addRunArtifact(tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&tests_run_cmd.step);
}

const playdate_target = blk: {
    @setEvalBranchQuota(6000);
    break :blk std.zig.CrossTarget.parse(.{
        .arch_os_abi = "thumb-freestanding-eabihf",
        .cpu_features = "cortex_m7-fp64-fp_armv8d16-fpregs64-vfp2-vfp3d16-vfp4d16",
    }) catch @compileError("could not parse Playdate compilation target strings");
};

const PlaydateGameExeOptions = struct {
    optimize: std.builtin.OptimizeMode,
    assets_dir: ?[]const u8 = null,
    modules: []const std.Build.ModuleDependency,
};

const PlaydateGameExe = struct {
    exe: *std.Build.CompileStep,
    dir: std.Build.LazyPath,
};

fn addPlaydateGameExe(
    b: *std.Build,
    root_source: std.Build.LazyPath,
    opts: PlaydateGameExeOptions,
) !PlaydateGameExe {
    const write_files = b.addWriteFiles();
    const source_dir = write_files.getDirectorySource();
    write_files.step.name = "write source directory";

    const lib = b.addSharedLibrary(.{
        .name = "pdex",
        .root_source_file = root_source,
        .optimize = opts.optimize,
        .target = .{},
    });
    for (opts.modules) |mod| {
        lib.addModule(mod.name, mod.module);
    }
    _ = write_files.addCopyFile(lib.getEmittedBin(), "pdex" ++ switch (os_tag) {
        .windows => ".dll",
        .macos => ".dylib",
        .linux => ".so",
        else => @panic("TODO: unsupported OS"),
    });

    const exe = b.addExecutable(.{
        .name = "pdex.elf",
        .root_source_file = root_source,
        .target = playdate_target,
        .optimize = opts.optimize,
    });
    for (opts.modules) |mod| {
        exe.addModule(mod.name, mod.module);
    }
    exe.force_pic = true;
    exe.link_emit_relocs = true;
    exe.setLinkerScript(.{ .path = "platform/playdate/link_map.ld" });
    if (opts.optimize == .ReleaseFast) {
        exe.omit_frame_pointer = true;
    }
    _ = write_files.addCopyFile(exe.getEmittedBin(), "pdex.elf");

    if (opts.assets_dir) |dir| {
        // WriteFile doesn't support copying whole directories.
        var assets = try b.build_root.handle.openIterableDir(dir, .{});
        defer assets.close();
        var iter = assets.iterate();
        while (try iter.next()) |entry| {
            if (entry.kind != .file) continue;
            const file_source = .{ .path = b.pathJoin(&.{ dir, entry.name }) };
            _ = write_files.addCopyFile(file_source, entry.name);
        }
    }

    return PlaydateGameExe{
        .exe = exe,
        .dir = source_dir,
    };
}

fn addPlaydateSimulatorRun(
    b: *std.Build,
    name: []const u8,
    game: PlaydateGameExe,
    playdate_sdk_path: ?[]const u8,
) !void {
    const default_sdk_path = switch (os_tag) {
        .windows => if (std.process.getEnvVarOwned(b.allocator, "USERPROFILE")) |home| path: {
            defer b.allocator.free(home);
            break :path std.fs.path.join(b.allocator, &.{ home, "Documents", "PlaydateSDK" }) catch @panic("OOM");
        } else |err| std.debug.panic("could not get %USERPROFILE%: {s}\n", .{@errorName(err)}),

        .macos => if (std.os.getenv("HOME")) |home|
            std.fs.path.join(b.allocator, &.{ home, "Developer", "PlaydateSDK" }) catch @panic("OOM")
        else
            @panic("could not get $HOME"),

        .linux => @panic("TODO: default Playdate SDK path on Linux"),
        else => @panic("TODO: unsupported OS"),
    };
    const sdk_path = playdate_sdk_path orelse default_sdk_path;

    const pdc_path = b.pathJoin(&.{ sdk_path, "bin", if (os_tag == .windows) "pdc.exe" else "pdc" });
    const pd_simulator_path = switch (os_tag) {
        .linux => b.pathJoin(&.{ sdk_path, "bin", "PlaydateSimulator" }),
        .macos => "open", // `open` focuses the window, while running the simulator directry doesn't.
        .windows => b.pathJoin(&.{ sdk_path, "bin", "PlaydateSimulator.exe" }),
        else => @panic("TODO: unsupported OS"),
    };

    const pdc = b.addSystemCommand(&.{ pdc_path, "--skip-unknown" });
    pdc.addDirectoryArg(game.dir);
    pdc.setName("pdc");
    const pdx_name = try std.mem.join(b.allocator, "", &.{ name, ".pdx" });
    const pdx = pdc.addOutputFileArg(pdx_name);

    b.installDirectory(.{
        .source_dir = pdx,
        .install_dir = .prefix,
        .install_subdir = pdx_name,
    });

    const run_playdate_simulator_cmd = b.addSystemCommand(&.{pd_simulator_path});
    run_playdate_simulator_cmd.addDirectoryArg(pdx);
    run_playdate_simulator_cmd.setName("PlaydateSimulator");

    const run_playdate_simulator_step = b.step(name, "Run example in the Playdate simulator");
    run_playdate_simulator_step.dependOn(&run_playdate_simulator_cmd.step);
    run_playdate_simulator_step.dependOn(b.getInstallStep());
}
