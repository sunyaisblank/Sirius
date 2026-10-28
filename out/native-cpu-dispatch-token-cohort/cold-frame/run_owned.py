"""Run the original frame controller only after current AT final admission."""
from pathlib import Path
from frame_authority import ao, load, check_acceptance

check_acceptance()
exec(compile(ao('run_owned.py'), str(Path.cwd() / 'out/native-cpu-v4-cohort/cold-frame/run_owned.py'), 'exec'))
