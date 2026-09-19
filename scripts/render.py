#!/usr/bin/env python3
"""Run explicitly selected example scenes; these outputs are not qualification."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import uuid


ROOT = Path(__file__).resolve().parents[1]

# Each entry contains only scene parameters. Workload, backend and output format
# are selected once on the command line, including for the former reference case.
SCENES = {
    "schwarzschild": ("Schwarzschild Page-Thorne disk", [
        "--metric", "Schwarzschild", "--distance", "30", "--inclination", "80",
        "--fov", "60", "--temperature-model", "NovikovThorne", "--disk-temperature", "50000"]),
    "kerr": ("Kerr spin 0.9 Page-Thorne disk", [
        "--metric", "Kerr", "--spin", "0.9", "--distance", "30", "--inclination", "80",
        "--fov", "60", "--temperature-model", "NovikovThorne", "--disk-temperature", "50000"]),
    "kerr-extremal": ("Near-extremal Kerr disk", [
        "--metric", "Kerr", "--spin", "0.998", "--distance", "20", "--inclination", "75",
        "--fov", "45", "--temperature-model", "NovikovThorne", "--disk-temperature", "50000"]),
    "schwarzschild-newtonian": ("Bounded Shakura-Sunyaev comparison", [
        "--metric", "Schwarzschild", "--distance", "30", "--inclination", "80",
        "--fov", "60", "--temperature-model", "ShakuraSunyaev", "--disk-temperature", "6500"]),
    "kerr-edge-on": ("Kerr disk at 89 degrees", [
        "--metric", "Kerr", "--spin", "0.9", "--distance", "30", "--inclination", "89", "--fov", "60"]),
    "kerr-face-on": ("Kerr disk at 15 degrees", [
        "--metric", "Kerr", "--spin", "0.9", "--distance", "30", "--inclination", "15", "--fov", "60"]),
    "kerr-volume": ("Kerr grey volume with procedural density", [
        "--metric", "Kerr", "--spin", "0.9", "--distance", "25", "--inclination", "75",
        "--fov", "50", "--volumetric", "--h-over-r", "0.15", "--turbulence"]),
    "schwarzschild-volume": ("Schwarzschild thick grey volume", [
        "--metric", "Schwarzschild", "--distance", "30", "--inclination", "80",
        "--fov", "60", "--volumetric", "--h-over-r", "0.2", "--turbulence"]),
    "kerr-film": ("Kerr cinematic display treatment", [
        "--metric", "Kerr", "--spin", "0.9", "--distance", "25", "--inclination", "75",
        "--fov", "50", "--cinematic", "--film", "--film-preset", "Interstellar"]),
    "kerr-extremal-film": ("Near-extremal Kerr cinematic close-up", [
        "--metric", "Kerr", "--spin", "0.998", "--distance", "15", "--inclination", "70",
        "--fov", "35", "--cinematic", "--film", "--film-preset", "Interstellar"]),
    "wormhole": ("Morris-Thorne throat", [
        "--metric", "Morris-Thorne", "--throat-radius", "1", "--distance", "15",
        "--inclination", "80", "--fov", "60", "--no-disk"]),
    "warp": ("Alcubierre moving wall", [
        "--metric", "Alcubierre", "--warp-velocity", "0.5", "--bubble-radius", "1",
        "--distance", "20", "--inclination", "80", "--fov", "90", "--no-disk"]),
}


def bounded_integer(minimum: int, maximum: int):
    def parse(value: str) -> int:
        number = int(value)
        if not minimum <= number <= maximum:
            raise argparse.ArgumentTypeError(f"must be between {minimum} and {maximum}")
        return number
    return parse


def preset_binary(preset: str, config: str) -> Path:
    """Resolve the declared preset, never search other builds for an executable."""
    presets = json.loads((ROOT / "CMakePresets.json").read_text(encoding="utf-8"))
    definitions = {item["name"]: item for item in presets["configurePresets"]}
    if preset not in definitions or definitions[preset].get("hidden", False):
        raise ValueError(f"unknown or hidden configure preset: {preset}")

    def resolve(name: str) -> dict:
        item = definitions[name]
        parents = item.get("inherits", [])
        if isinstance(parents, str):
            parents = [parents]
        inherited = {}
        # CMake gives the first inherited preset precedence over later ones.
        for parent in reversed(parents):
            inherited.update(resolve(parent))
        inherited.update(item)
        return inherited

    selected = resolve(preset)
    directory = selected["binaryDir"].replace("${sourceDir}", str(ROOT))
    directory = directory.replace("${presetName}", preset)
    if "$" in directory:
        raise ValueError("unsupported preset binaryDir expansion; use --binary")
    binary = Path(directory) / "src/sirius/app"
    generator = selected["generator"]
    if generator.startswith("Visual Studio") or generator in {"Xcode", "Ninja Multi-Config"}:
        binary /= config
    windows = generator.startswith("Visual Studio") or os.name == "nt"
    return binary / ("sirius.exe" if windows else "sirius")


def plan(args: argparse.Namespace, run_id: str) -> dict:
    binary = args.binary.resolve() if args.binary else preset_binary(args.preset, args.config)
    image_dir = ROOT / "renders" / run_id
    log_dir = ROOT / "out" / "render" / run_id
    cases = list(SCENES) if args.all else list(dict.fromkeys(args.scene))
    commands = []
    for name in cases:
        output = image_dir / f"{name}.{args.format}"
        command = [str(binary), "render", *SCENES[name][1],
                   "--width", str(args.width), "--height", str(args.height),
                   "--samples", str(args.samples), "--backend", args.backend,
                   "--output", str(output)]
        commands.append({"scene": name, "command": command, "output": str(output),
                         "log": str(log_dir / f"{name}.log")})
    return {"kind": "sirius-example-render", "binary": str(binary),
            "image_directory": str(image_dir), "log_directory": str(log_dir),
            "status": "planned", "cases": commands}


def execute(document: dict) -> int:
    binary = Path(document["binary"])
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError(f"Sirius binary is not executable: {binary}; build the selected preset first")
    image_dir = Path(document["image_directory"])
    log_dir = Path(document["log_directory"])
    # A new directory is mandatory: an old image can never satisfy a new run.
    image_dir.mkdir(parents=True, exist_ok=False)
    try:
        log_dir.mkdir(parents=True, exist_ok=False)
    except OSError:
        image_dir.rmdir()
        raise
    manifest = log_dir / "run.json"

    def record():
        pending = manifest.with_suffix(".tmp")
        pending.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
        pending.replace(manifest)

    document["status"] = "running"
    record()
    try:
        for case in document["cases"]:
            print(f"Rendering {case['scene']} (log: {case['log']})", flush=True)
            case["status"] = "running"
            record()
            with Path(case["log"]).open("w", encoding="utf-8") as log:
                result = subprocess.run(case["command"], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            case["returncode"] = result.returncode
            output = Path(case["output"])
            if result.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
                case["status"] = "failed"
                output.unlink(missing_ok=True)
                document["status"] = "failed"
                print(f"Render failed or emitted no image; see {case['log']}", file=sys.stderr)
                return 1
            case["status"] = "complete"
            record()
        document["status"] = "complete"
        print(f"Images: {image_dir}\nLogs and command record: {log_dir}")
        return 0
    except (OSError, KeyboardInterrupt):
        document["status"] = "interrupted"
        for case in document["cases"]:
            if case.get("status") == "running":
                case["status"] = "interrupted"
                Path(case["output"]).unlink(missing_ok=True)
        raise
    finally:
        record()
        if not any(image_dir.iterdir()):
            image_dir.rmdir()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--scene", action="append", choices=SCENES, default=[],
                           help="scene to render; repeat for a selected batch")
    selection.add_argument("--all", action="store_true", help="explicitly render every example")
    default_preset = {"win32": "windows-msvc", "darwin": "macos"}.get(sys.platform, "linux-gcc")
    parser.add_argument("--preset", default=default_preset,
                        help=f"exact configure preset (default: {default_preset})")
    parser.add_argument("--config", choices=("Debug", "Release", "RelWithDebInfo", "MinSizeRel"),
                        default="Release", help="configuration for multi-config presets")
    parser.add_argument("--binary", type=Path, help="explicit executable instead of the preset binary")
    parser.add_argument("--width", type=bounded_integer(128, 8192), default=128)
    parser.add_argument("--height", type=bounded_integer(128, 8192), default=128)
    parser.add_argument("--samples", type=bounded_integer(1, 4096), default=1)
    parser.add_argument("--backend", choices=("cpu", "vulkan", "auto"), default="cpu")
    parser.add_argument("--format", choices=("png", "ppm", "exr"), default="png")
    parser.add_argument("--dry-run", action="store_true", help="print exact commands without creating files")
    args = parser.parse_args()
    if not args.scene and not args.all:
        for name, (description, _) in SCENES.items():
            print(f"{name:26} {description}")
        print("Select --scene NAME (repeatable) or --all to render; --dry-run previews commands.")
        return 0
    try:
        run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
        document = plan(args, run_id)
        if args.dry_run:
            print(json.dumps(document, indent=2))
            return 0
        return execute(document)
    except (OSError, ValueError) as error:
        print(f"render workflow: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Render interrupted; partial image removed and logs retained.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
