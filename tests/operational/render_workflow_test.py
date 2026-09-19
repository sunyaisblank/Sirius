#!/usr/bin/env python3
"""Exercise the operator workflow with child processes that never render."""

import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "render.py"
SPEC = importlib.util.spec_from_file_location("render_workflow", SCRIPT)
workflow = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(workflow)
RUN = subprocess.run


class RenderWorkflowTest(unittest.TestCase):
    def setUp(self):
        (ROOT / "out").mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix="render-workflow-", dir=ROOT / "out")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "checkout with spaces"
        self.root.mkdir()
        (self.root / "CMakePresets.json").write_bytes((ROOT / "CMakePresets.json").read_bytes())
        self.addCleanup(patch.stopall)
        patch.object(workflow, "ROOT", self.root).start()
        self.args = argparse.Namespace(binary=Path(sys.executable), preset="linux-gcc",
                                      config="Release", scene=["kerr", "wormhole"], all=False,
                                      width=128, height=128, samples=1, backend="cpu", format="png")
        self.output = io.StringIO()
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.stack.enter_context(contextlib.redirect_stdout(self.output))
        self.stack.enter_context(contextlib.redirect_stderr(self.output))

    def run_fake(self, body):
        fake = self.root / "fake renderer.py"
        fake.write_text("import pathlib, sys\n" + body, encoding="utf-8")
        document = workflow.plan(self.args, "example")
        # Real child exit codes and bytes exercise the workflow without needing
        # an executable shell script (which would exclude native Windows).
        def child(command, **kwargs):
            return RUN([sys.executable, str(fake), *command[1:]], **kwargs)
        with patch.object(workflow.subprocess, "run", side_effect=child) as calls:
            result = workflow.execute(document)
        manifest = json.loads((self.root / "out/render/example/run.json").read_text(encoding="utf-8"))
        return result, manifest, calls.call_count

    def test_success_records_commands_and_keeps_logs_out_of_images(self):
        result, manifest, count = self.run_fake(
            "pathlib.Path(sys.argv[sys.argv.index('--output') + 1]).write_bytes(b'image')\n"
            "print('renderer log')\n")
        self.assertEqual((result, manifest["status"], count), (0, "complete", 2))
        self.assertEqual(sorted(p.name for p in (self.root / "renders/example").iterdir()),
                         ["kerr.png", "wormhole.png"])
        for case in manifest["cases"]:
            self.assertEqual(case["returncode"], 0)
            self.assertEqual(Path(case["log"]).read_text().strip(), "renderer log")
            self.assertIn(str(Path(case["output"])), case["command"])

    def test_renderer_failure_removes_partial_output_and_stops_batch(self):
        result, manifest, count = self.run_fake(
            "pathlib.Path(sys.argv[sys.argv.index('--output') + 1]).write_bytes(b'partial')\n"
            "print('failure detail')\n"
            "sys.exit(19)\n")
        self.assertEqual((result, manifest["status"], count), (1, "failed", 1))
        self.assertEqual(manifest["cases"][0]["returncode"], 19)
        self.assertFalse((self.root / "renders/example").exists())
        self.assertIn("failure detail", Path(manifest["cases"][0]["log"]).read_text())

    def test_exit_zero_without_an_image_is_failure(self):
        result, manifest, count = self.run_fake("sys.exit(0)\n")
        self.assertEqual((result, manifest["status"], count), (1, "failed", 1))

    def test_exit_zero_with_an_empty_image_is_failure(self):
        result, manifest, count = self.run_fake(
            "pathlib.Path(sys.argv[sys.argv.index('--output') + 1]).touch()\n")
        self.assertEqual((result, manifest["status"], count), (1, "failed", 1))

    def test_existing_images_cannot_satisfy_a_new_run(self):
        document = workflow.plan(self.args, "example")
        output = Path(document["cases"][0]["output"])
        output.parent.mkdir(parents=True)
        output.write_bytes(b'old image')
        with patch.object(workflow.subprocess, "run") as child:
            with self.assertRaises(FileExistsError):
                workflow.execute(document)
            child.assert_not_called()
        self.assertEqual(output.read_bytes(), b'old image')

    def test_interruption_removes_partial_image_but_retains_record(self):
        document = workflow.plan(self.args, "example")
        def interrupt(*args, **kwargs):
            Path(document["cases"][0]["output"]).write_bytes(b'partial')
            raise KeyboardInterrupt
        with patch.object(workflow.subprocess, "run", side_effect=interrupt):
            with self.assertRaises(KeyboardInterrupt):
                workflow.execute(document)
        manifest = json.loads((self.root / "out/render/example/run.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["status"], "interrupted")
        self.assertFalse(Path(document["cases"][0]["output"]).exists())

    def test_missing_binary_is_reported_without_creating_outputs(self):
        self.args.binary = self.root / "missing"
        with self.assertRaisesRegex(ValueError, "not executable"):
            workflow.execute(workflow.plan(self.args, "example"))
        self.assertFalse((self.root / "renders").exists())
        self.assertFalse((self.root / "out").exists())

    def test_preset_resolution_uses_exact_configuration(self):
        name = "sirius.exe" if sys.platform == "win32" else "sirius"
        self.assertEqual(workflow.preset_binary("linux-gcc", "Debug"),
                         self.root / "bin/linux-gcc/src/sirius/app" / name)
        self.assertEqual(workflow.preset_binary("windows-msvc", "Debug"),
                         self.root / "bin/windows-msvc/src/sirius/app/Debug/sirius.exe")
        for preset in ("base", "unknown", "../linux-gcc"):
            with self.assertRaises(ValueError):
                workflow.preset_binary(preset, "Release")

    def test_selection_and_dry_run_never_start_the_renderer(self):
        with patch.object(sys, "argv", [str(SCRIPT)]), patch.object(workflow, "execute") as execute:
            self.assertEqual(workflow.main(), 0)
            execute.assert_not_called()
        self.assertIn("Select --scene", self.output.getvalue())
        self.output.seek(0)
        self.output.truncate()
        with patch.object(sys, "argv", [str(SCRIPT), "--all", "--dry-run"]), \
                patch.object(workflow, "execute") as execute:
            self.assertEqual(workflow.main(), 0)
            execute.assert_not_called()
        document = json.loads(self.output.getvalue())
        self.assertEqual(len(document["cases"]), 12)
        self.assertFalse((self.root / "renders").exists())
        self.assertFalse((self.root / "out").exists())


if __name__ == "__main__":
    unittest.main()
