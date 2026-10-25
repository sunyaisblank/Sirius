"""Apply the established installed-zlib transition to this fresh provider cohort."""
from pathlib import Path
basis = Path.cwd() / 'out/native-cpu-zero-term-cohort/runtime_transition.py'
source = basis.read_text().replace('out/native-cpu-zero-term/provider', 'out/native-cpu-field-simd/provider')
source = source.replace('the S provider', 'the T provider').replace('the S provider', 'the T provider')
exec(compile(source, str(basis), 'exec'))
