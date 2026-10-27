"""Keep the original frame verdict, XML, resources and terminal-owner reader."""
from pathlib import Path
from frame_authority import ao, load, check_acceptance

check_acceptance()
exec(compile(ao('read_result.py'), str(Path.cwd() / 'out/native-cpu-v4-cohort/cold-frame/read_result.py'), 'exec'))
