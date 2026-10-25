"""Reuse the original one-object strict adapter compiler runner."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-simd-rows/build_adapter.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '3dc7d9194c09171ef4d1f4b9be8ff4bd0866d07a68baff09f2a401ad1bcf4236'
exec(compile(basis.read_text(), str(basis), 'exec'))
