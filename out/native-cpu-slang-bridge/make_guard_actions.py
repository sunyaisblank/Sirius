"""Seal three finite original4B fallback captures; never execute them here."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
O = W / 'ordinary'
B = R / 'bin/linux-gcc/native-retained-ownership'


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


assert not (O / 'actions.json').exists()
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
build = json.loads((W / 'provider/provider-build-result.json').read_text())
assert json.loads((W / 'provider/build/owner.json').read_text())['passed']
assert json.loads((W / 'provider/linked-native-proof.json').read_text())['passed']
library, icd = R / build['library'], R / build['ICD']
assert pin(library) == build['library_identity'] and pin(icd) == build['ICD_identity']
inputs = {str(path.relative_to(R)) for path in O.iterdir() if path.is_file()}
inputs |= {str(path.relative_to(R)) for path in B.rglob('*') if path.is_file()}
inputs |= {str(path.relative_to(R)) for path in
           (R / 'out/native-original-input-ownership/reused/prior-proof').rglob('*') if path.is_file()}
inputs |= {str(path.relative_to(R)) for path in (W / 'provider').iterdir() if path.is_file()}
inputs |= {'out/lavapipe-acqrel-ordering/ordinary/capture.py',
           'out/native-cpu-row-ownership/cold-frame/run_owned.py',
           'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
           'out/native-cpu-row-ownership/compiler-result.json',
           str(library.relative_to(R)), str(icd.relative_to(R)),
           str(W / 'native_observer.py'), str(W / 'prepare_guards.py'),
           str(W / 'make_guard_actions.py'), str(W / 'read_linked_native.py')}
for item in build['native_objects']:
    inputs.add(item['path'])
modules = json.loads((R / 'out/native-cpu-row-ownership/compiler-result.json').read_text())['modules']
for stage in ('Transport', 'Endpoint', 'DopriPhase'):
    item = modules[stage]
    module = R / item['module']
    assert pin(module) == {key: item[key] for key in ('bytes', 'sha256')}
    inputs |= {str(module.relative_to(R)), str(module.with_suffix('.spvasm').relative_to(R))}
for path in [library, B / 'short_input_jit_probe', B / 'optimizer_observer.so',
             Path('/usr/bin/gdb'), Path('/usr/bin/python3'), Path('/usr/bin/nm'),
             Path('/usr/bin/objdump'), Path('/bin/bash')]:
    inputs |= {str(path), str(path.resolve(strict=True))}
    dependencies = subprocess.check_output(['ldd', str(path)], text=True)
    for name in re.findall(r'(/[^\s()]+)', dependencies):
        dependency = Path(name)
        if dependency.is_file():
            inputs |= {str(dependency), str(dependency.resolve(strict=True))}
inputs |= {'/etc/ld.so.cache', '/usr/bin/ldd', str(Path('/usr/bin/ldd').resolve(strict=True))}
assert all((R / name).is_file() for name in inputs)
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
actions = {}
for slug, stage in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
    actions[slug] = {
        'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
        'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'argv': ['/usr/bin/python3', '-B', str(O / 'capture.py'), stage, '{folder}'],
        'frozen_inputs': sorted(inputs), 'software_provider': True,
        'runtime_icd': str(icd), 'runtime_provider': str(library),
        'scope': 'One original4B guard/canary dispatch per exact module: passive native entry+returnedfalse/outputunchanged then real ordinary JIT entry/features/IEEE32 callback restoration. Static native object proof already passed. Negative fallback/code boundary only; no successful native arithmetic, numerical/benefit/frame/release claim.',
    }
(O / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'sealed_cases': len(actions), 'sealed_inputs': len(inputs)}))
