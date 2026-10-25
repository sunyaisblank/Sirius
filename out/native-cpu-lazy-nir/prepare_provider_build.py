"""Reuse the existing incremental planner and original resource limits."""
from pathlib import Path
basis = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_build.py'
exec(compile(basis.read_text(), str(basis), 'exec'))
