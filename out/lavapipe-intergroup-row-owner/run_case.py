"""Reuse the original acceptance authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-acqrel-ordering/run_case.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'f9c87d6108ffd8ce17d2cf9867e0ab6340baef65310b6437eb6407a321c5b966'
_source = _basis.read_text()
assert 'bin/linux-gcc/native-original-input-ownership/matched/audited_backend_tests' in _source
_source = _source.replace('bin/linux-gcc/native-original-input-ownership/matched/audited_backend_tests', 'bin/linux-gcc/native-cpu-row-ownership/matched/audited_backend_tests')
exec(compile(_source, str(_basis), "exec"))
