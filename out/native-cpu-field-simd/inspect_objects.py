"""Use the existing ISA reader for the new layer exports and actual receipts."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-zero-term/inspect_objects.py'
source = basis.read_text()
assert hashlib.sha256(source.encode()).hexdigest() == 'f7153a08cab7ef084e79e5b63392934803e47d690cb31566e1caa918b0ee5288'
before = 'exec(compile(source, str(basis), \'exec\'))'
assert source.count(before) == 1
source = source.replace(before, '''source = source.replace('ispc-zero-compile-result.json', 'ispc-compile-result.json')
source = source.replace('compile/ispc-zero/owner.json', 'compile/ispc/owner.json')
source = source.replace("+ '_packet'", "+ '_field_layer'")
source = source.replace('object-code-inspection-zero.json', 'object-code-inspection.json')
exec(compile(source, str(basis), 'exec'))''')
exec(compile(source, str(basis), 'exec'))
