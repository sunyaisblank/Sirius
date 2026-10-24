"""Compile the sealed adapter once; keep all original object bytes unchanged."""
from pathlib import Path
import hashlib
import json
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
prepared = json.loads((work / 'prepared-adapter.json').read_text())


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


for name, value in prepared['inputs'].items():
    assert pin(root / name) == value, name
output = root / prepared['output']
assert not output.exists() and not (work / 'adapter-compile-result.json').exists()
result = subprocess.run(prepared['argv'])
record = {'objects': [{'exit': result.returncode, 'argv': prepared['argv']}],
          'source_inputs': prepared['inputs'], 'all_prepared_source_bytes_unchanged': True}
for name, value in prepared['inputs'].items():
    assert pin(root / name) == value, name
if output.is_file():
    record['objects'][0]['object'] = {'path': prepared['output'], **pin(output)}
(work / 'adapter-compile-result.json').write_text(json.dumps(record, indent=2) + '\n')
raise SystemExit(result.returncode)
