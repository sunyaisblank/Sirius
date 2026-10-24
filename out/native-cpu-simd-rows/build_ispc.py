"""Compile sealed ISPC sources serially once; preserve all diagnostics."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
prepared = json.loads((work / 'prepared-ispc-uniform-inputs.json').read_text())
receipt = work / 'ispc-uniform-compile-result.json'
destination = root / 'bin/linux-gcc/native-cpu-simd-rows/objects-uniform'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def verify():
    assert not os.environ.get('ISPC_ARGS', '')
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == prepared['revision']
    for name, identity in prepared['inputs'].items():
        assert pin(root / name) == identity, name


def save():
    receipt.write_text(json.dumps(result, indent=2) + '\n')


verify()
assert not receipt.exists() and not destination.exists()
destination.mkdir(parents=True)
result = {'revision': prepared['revision'], 'objects': [], 'all_objects_compiled': False,
          'source_inputs': prepared['inputs'], 'scope': prepared['scope']}
with receipt.open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
exit_code = 0
for item in prepared['objects']:
    assert not (root / item['output']).exists()
    row = {'stage': item['stage'], 'argv': item['argv'], 'status': 'started'}
    result['objects'].append(row)
    save()
    print(json.dumps(row), flush=True)
    process = subprocess.run(item['argv'], cwd=root)
    row.update(status='terminal', exit=process.returncode)
    output = root / item['output']
    if output.is_file():
        row['object'] = {'path': item['output'], **pin(output)}
    save()
    if process.returncode:
        exit_code = process.returncode
        break
verify()
result['all_prepared_source_bytes_unchanged'] = True
result['all_objects_compiled'] = len(result['objects']) == 3 and all(
    row['exit'] == 0 and row.get('object', {}).get('bytes', 0) > 0 for row in result['objects'])
save()
raise SystemExit(exit_code or (0 if result['all_objects_compiled'] else 1))
