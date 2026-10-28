"""Reuse the original owned-session, memory and cleanup guard unchanged."""
from pathlib import Path
from frame_authority import ao, load

exec(compile(ao('owned_guard.py'), str(Path.cwd() / 'out/native-cpu-v4-cohort/cold-frame/owned_guard.py'), 'exec'))
