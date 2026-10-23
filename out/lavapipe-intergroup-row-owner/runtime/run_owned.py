"""Reuse the original owner and limits with this runtime action folder."""
from pathlib import Path
_preserved = Path(__file__).resolve().parents[2] / 'native-cpu-row-ownership/cold-frame/run_owned.py'
exec(compile(_preserved.read_text(), str(_preserved), 'exec'))
