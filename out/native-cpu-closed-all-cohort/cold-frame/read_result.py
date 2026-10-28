"""Read the sole terminal frame, including failure; never retry it here."""
from pathlib import Path
from frame_authority import ao, load, check_acceptance

check_acceptance()
exec(compile(ao('read_result.py'), str(Path.cwd() / 'out/native-cpu-v4-cohort/cold-frame/read_result.py'), 'exec'))
