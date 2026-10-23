"""Reconfigure the existing private Mesa build for the three sealed native objects."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
B = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
P = W / 'provider'


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


contract = json.loads((P / 'source-contract.json').read_text())
review = json.loads((P / 'review-admission.json').read_text())
assert review['passed'] and review['source_contract'] == pin(P / 'source-contract.json')
assert review['patch'] == pin(P / 'provider.patch')
source = {str(path.relative_to(S)): pin(path) for path in S.rglob('*') if path.is_file()}
encoded = json.dumps(source, sort_keys=True, separators=(',', ':')).encode()
assert {'bytes': len(encoded), 'sha256': hashlib.sha256(encoded).hexdigest()} == contract['full_candidate_source_identity']
for item in contract['native_objects']:
    assert pin(R / item['path']) == {key: item[key] for key in ('bytes', 'sha256')}
argv = [str(R / 'bin/toolchains/python/bin/meson'), 'setup', '--reconfigure', str(B), str(S)]
process = subprocess.run(argv)
(P / 'configure-result.json').write_text(json.dumps({
    'argv': argv, 'exit': process.returncode,
    'existing_build_directory': str(B.relative_to(R)),
    'options_changed': False, 'new_toolchain': False,
    'scope': 'Existing private Mesa configuration regenerated to link3 already-sealed native objects. No Ninja build/provider loading/numerical/performance result.',
}, indent=2) + '\n')
raise SystemExit(process.returncode)
