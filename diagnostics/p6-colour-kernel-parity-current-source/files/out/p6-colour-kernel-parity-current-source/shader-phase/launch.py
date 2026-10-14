"""Root-only sequential launch through the unchanged prior process owner."""
from pathlib import Path
import os
import sys

sys.dont_write_bytecode = True
from check import reader

phase = Path(__file__).resolve().parent
assert len(sys.argv) == 2 and sys.argv[1] in {"cie", "xyz"}
case = sys.argv[1]
original, _, _ = reader()
original["preflight"](case)
assert not (phase / case).exists(), "No retry or overwrite"
assert not (phase / ("accepted-" + case + ".json")).exists()
if case == "xyz":
    original["validate"]("cie")
os.execv(sys.executable, [sys.executable, str(phase / "run_owned.py"), case])
