"""Plan only dirty C outputs; reuse X's unchanged configuration and native objects."""
from pathlib import Path
import hashlib
import json
import re
import shlex
import subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
B = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
X = R / 'out/native-cpu-predicated-arithmetic/provider'
doc = lambda p: json.loads(p.read_text())
pin = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert doc(P / 'review-admission.json')['passed']
assert doc(X / 'configure/owner.json')['passed'] and doc(X / 'build/owner.json')['passed']
contract = doc(P / 'source-contract.json')
assert contract['new_native_compilations'] == 0
for name, expected in contract['source_inputs'].items():
    assert pin(R / name) == expected, name
source = {str(p.relative_to(S)): pin(p) for p in S.rglob('*') if p.is_file()}
encoded = json.dumps(source, sort_keys=True, separators=(',', ':')).encode()
assert {'bytes': len(encoded), 'sha256': hashlib.sha256(encoded).hexdigest()} == contract['full_candidate_source_identity']
dry = subprocess.check_output(['ninja', '-C', str(B), '-n', '-v', 'src/gallium/targets/lavapipe/libvulkan_lvp.so'], text=True)
commands = [m[1] for line in dry.splitlines() if (m := re.match(r'^\[\d+/\d+\] (.*)$', line))]
assert commands
mutable = set()
for command in commands:
    words = shlex.split(command)
    outputs = {words[i + 1] for i, word in enumerate(words[:-1]) if word in {'-o', '--output', 'csr', 'csrD', 'csrDT', 'crs', 'crsD', 'rcs', 'rcsD'}}
    assert outputs, ('Review unsupported build command', command)
    for output in outputs:
        path = (B / output).resolve(); assert path.is_relative_to(B)
        mutable.add(str(path.relative_to(R)))
target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'
assert str(target.relative_to(R)) in mutable
reused = {str(p.relative_to(R)): pin(p) for p in B.rglob('*') if p.is_file()
          and str(p.relative_to(R)) not in mutable and p.name not in {'.ninja_log', '.ninja_deps', '.ninja_lock'}}
plan = {'dry_run_commands': commands, 'mutable_outputs': sorted(mutable), 'reused_build_inputs': reused,
    'previous_target': pin(target), 'configuration_reused_without_reconfigure': {
        'owner': pin(X / 'configure/owner.json'), 'result': pin(X / 'configure-result.json'),
        'build_ninja': pin(B / 'build.ninja'), 'reason': 'Only four C source files change. Same native object paths/Meson files/options/toolchain.'}}
assert not (P / 'incremental-plan.json').exists() and not (P / 'actions.json').exists()
(P / 'incremental-plan.json').write_text(json.dumps(plan, indent=2) + '\n')
(P / 'incremental-dry-run.txt').write_text(dry)
base = doc(X / 'actions.json')['build']
inputs = {name for name in base['frozen_inputs'] if not (R / name).is_relative_to(B)}
inputs |= set(contract['source_inputs']) | set(reused)
inputs |= {str(p.relative_to(R)) for p in S.rglob('*') if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()}
inputs |= {row['path'] for row in contract['native_objects']}
inputs |= {str(p.relative_to(R)) for p in (X / 'configure/owner.json', X / 'configure-result.json', X / 'build/owner.json')}
assert all((R / name).is_file() for name in inputs)
action = {**base, 'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(W / 'build_provider.py')], 'frozen_inputs': sorted(inputs),
    'scope': 'One incremental Y provider build: four changed C sources, exact X configuration/build cache/native8 reused. Zero Meson reconfiguration, native/host/fixture compilation or new toolchain. Same limits and strict options; no runtime/performance/scientific/release qualification.'}
(P / 'actions.json').write_text(json.dumps({'build': action}, indent=2) + '\n')
for name in ('run_owned.py', 'owned_guard.py'):
    (P / name).write_bytes((W / name).read_bytes())
print(json.dumps({'pending_commands': len(commands), 'reused_build_files': len(reused), 'reconfigure': False, 'new_native_compilations': 0}))
