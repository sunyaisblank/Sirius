"""Compatibility entry point for the retained camera witness generator."""
import sys
sys.dont_write_bytecode = True

import importlib.util
from pathlib import Path

_source = Path(__file__).resolve().parents[3] / "src/sirius/kernels/retained_program.py"
_spec = importlib.util.spec_from_file_location("sirius_retained_program", _source)
_program = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_program)
build_program = _program.build_camera_program
