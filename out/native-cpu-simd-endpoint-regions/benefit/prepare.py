"""Reuse original ABBA cases, counts, limits and numerical/Decimal gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'ba3c35750f9e1b0b8e876b084d1f3f738ab34d1080709080dd5d4dc76915a440'
source = basis.read_text()
source = source.replace('current B', 'current L').replace('native B', 'native L')
exec(compile(source, str(basis), 'exec'))
