"""Compile the three sealed typed adapters once; never execute native shader code."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
prepared = json.loads((W / 'prepared-batch-inputs.json').read_text())


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def verify():
    for name, identity in prepared['inputs'].items():
        assert pin(R / name) == identity, name


verify()
results = []
exit_code = 0
for item in prepared['objects']:
    output = R / item['output']
    assert not output.exists()
    process = subprocess.run(item['argv'])
    result = {'stage': item['stage'], 'argv': item['argv'], 'exit': process.returncode}
    if output.exists():
        result['object'] = {'path': item['output'], **pin(output)}
    results.append(result)
    if process.returncode:
        exit_code = process.returncode
        break
verify()
(W / 'batch-compile-result.json').write_text(json.dumps({
    'objects': results, 'all_prepared_source_bytes_unchanged': True,
    'source_inputs': prepared['inputs'],
    'scope': 'Strict native C++26 compile-only typedABI/packet adapters. No object execution, driver integration, numerical/FP equivalence or benefit qualification.',
}, indent=2) + '\n')
raise SystemExit(exit_code)
