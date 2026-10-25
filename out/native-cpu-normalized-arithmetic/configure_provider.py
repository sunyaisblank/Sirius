"""Reuse the unchanged private Mesa configuration command and source checks."""
from pathlib import Path
basis = Path.cwd() / 'out/native-cpu-slang-bridge/configure_provider.py'
exec(compile(basis.read_text(), str(basis), 'exec'))
