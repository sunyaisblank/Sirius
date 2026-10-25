"""Reuse original ABBA cases, counts, limits and numerical/Decimal gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/readers/decimal_check.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '01a969846717bd6ff6d4a5637cf8548801426d79a14145a3d305109885e4997d'
source = basis.read_text()
assert 'out/native-cpu-active-prefix' in source
source = source.replace('out/native-cpu-active-prefix', 'out/native-cpu-task-cardinality-cohort')
exec(compile(source, str(basis), 'exec'))
