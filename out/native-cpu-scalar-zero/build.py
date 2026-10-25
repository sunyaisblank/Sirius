"""Compile the selected sealed three objects once, retaining partial failures."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
assert sys.argv[1:] in (['cpp'], ['ispc'])
kind = sys.argv[1]
prepared = json.loads((W / 'prepared-compile.json').read_text())
receipt = W / (kind + '-compile-result.json')


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def verify():
    assert not os.environ.get('ISPC_ARGS')
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == prepared['revision']
    for name, expected in prepared['inputs'].items():
        assert pin(R / name) == expected, name


def save():
    receipt.write_text(json.dumps(result, indent=2) + '\n')


verify()
assert not receipt.exists()
items = [item for item in prepared['objects'] if item['kind'] == kind]
assert len(items) == 3 and all(not (R / item['output']).exists() for item in items)
(R / 'bin/linux-gcc/native-cpu-scalar-zero/objects').mkdir(parents=True, exist_ok=True)
result = {'revision': prepared['revision'], 'kind': kind, 'objects': [],
          'source_inputs': prepared['inputs'], 'scope': prepared['scope'], 'all_objects_compiled': False}
with receipt.open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
code = 0
for item in items:
    row = {'stage': item['stage'], 'argv': item['argv'], 'status': 'started'}
    result['objects'].append(row); save()
    print(json.dumps(row), flush=True)
    code = subprocess.run(item['argv'], cwd=R).returncode
    row.update(argv_executed=True, exit=code, status='terminal')
    output = R / item['output']
    if output.is_file():
        row['object'] = {'path': item['output'], **pin(output)}
    save()
    if code:
        break
verify()
result['new_compilations'] = len(result['objects'])
result['all_prepared_source_bytes_unchanged'] = True
result['all_objects_compiled'] = len(result['objects']) == 3 and all(
    item['exit'] == 0 and item.get('object', {}).get('bytes', 0) > 0 for item in result['objects'])
save()
raise SystemExit(code or (0 if result['all_objects_compiled'] else 1))
