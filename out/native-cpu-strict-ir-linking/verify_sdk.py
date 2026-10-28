"""Run the unchanged extraction/version logic after the completed download."""
from pathlib import Path
import hashlib
import textwrap

root = Path.cwd()
basis = root / 'out/native-cpu-strict-ir-linking/acquire.py'
source = basis.read_text()
assert hashlib.sha256(source.encode()).hexdigest() == '9533b67e4c34a61a5fe38a12eab058cb200e343e86ad68bf0f7abf2d2a3a732c'
scope = {'__file__': str(basis), '__name__': 'verified_sdk_extraction'}
prefix, = source.split('phase, = sys.argv[1:]')[:1]
exec(compile(prefix, str(basis), 'exec'), scope)
archive = scope['B'] / scope['ASSET']
assert scope['pin'](archive)['sha256'] == scope['SHA']
assert not scope['D'].exists()
scope['archive'] = archive
start = source.index('    D.mkdir()\n')
end = source.rindex('\nelse:\n')
exec(compile(textwrap.dedent(source[start:end]), str(basis), 'exec'), scope)
