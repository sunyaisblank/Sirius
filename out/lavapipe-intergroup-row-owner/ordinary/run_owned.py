"""Use the original process owner with this action folder."""
from pathlib import Path
P = Path(__file__).resolve().parents[2] / 'native-cpu-row-ownership/cold-frame/run_owned.py'
exec(compile(P.read_text(), str(P), 'exec'))
