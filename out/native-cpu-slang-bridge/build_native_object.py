"""Compile exactly the already sealed source/header pair; never run the object."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
prepared = json.loads((W / 'prepared-native-inputs.json').read_text())
def pin(p):
    data = p.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def verify():
    for name, identity in prepared['inputs'].items():
        assert pin(R / name) == identity, name
verify()
output = R / prepared['output']
assert not output.exists()
result = subprocess.run(prepared['argv'])
verify()
record = {'argv': prepared['argv'], 'exit': result.returncode,
          'all_prepared_source_bytes_unchanged': True,
          'source_inputs': prepared['inputs'],
          'scope': 'One strict-policy native C++26 object/typed-ABI screen; no object execution, provider integration, FP/runtime equivalence or benefit qualification.'}
if output.exists():
    record['actual_object'] = {'path': str(output.relative_to(R)), **pin(output)}
(W / 'compile-result.json').write_text(json.dumps(record, indent=2) + '\n')
raise SystemExit(result.returncode)
