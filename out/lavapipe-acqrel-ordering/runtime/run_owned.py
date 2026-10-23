"""Reuse the preserved owner with this action directory."""
from pathlib import Path
_preserved_runner = Path(__file__).resolve().parents[2] / "native-cpu-row-ownership/cold-frame/run_owned.py"
exec(compile(_preserved_runner.read_text(), str(_preserved_runner), "exec"))
