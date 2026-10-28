"""Original H frame/resources only after current AX all8/all22/Decimal."""
from pathlib import Path
from frame_authority import ao, load, check_acceptance

check_acceptance()
exec(compile(ao('prepare.py'), str(Path.cwd() / 'out/native-cpu-v4-cohort/cold-frame/prepare.py'), 'exec'))
