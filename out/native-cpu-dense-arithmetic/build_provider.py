"""Reuse the current incremental Mesa build; stage only the new private provider."""
from pathlib import Path
basis = Path.cwd() / 'out/native-cpu-slang-bridge/build_provider.py'
source = basis.read_text()
before = 'native-cpu-slang-bridge/provider'
assert source.count(before) == 1
source = source.replace(before, 'native-cpu-dense-arithmetic/provider')
exec(compile(source, str(basis), 'exec'))
