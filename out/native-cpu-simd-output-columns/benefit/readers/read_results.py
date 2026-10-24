"""Reuse original ABBA cases, counts, limits and numerical/Decimal gates."""
from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-active-prefix/benefit/readers/read_results.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '1371748434db9de22af3555de6e2b27b25e9bd25f7928d1937f8f0202a7ace92'
source = basis.read_text()
assert 'out/native-cpu-active-prefix' in source
source = source.replace('out/native-cpu-active-prefix', 'out/native-cpu-simd-output-columns')
exec(compile(source, str(basis), 'exec'))
