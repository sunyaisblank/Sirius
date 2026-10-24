"""Reuse the exact existing private Mesa configuration and options."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-slang-bridge/configure_provider.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '887a9d6f484f9791aea3ba297899b3275b84ceca726fdfd6dc83d296fef30c85'
exec(compile(basis.read_text(), str(basis), 'exec'))
