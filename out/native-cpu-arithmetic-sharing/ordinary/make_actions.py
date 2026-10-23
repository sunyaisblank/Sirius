"""Seal the two original ordinary-code guard captures; never run them here."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W.parent
B = R / 'bin/linux-gcc/native-retained-ownership'
def pin(p):
    data = p.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def key(p):
    return str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
assert not (W / 'actions.json').exists()
assert not any((W / name).exists() for name in ['endpoint', 'dopri_phase'])
result = json.loads((P / 'provider/provider-build-result.json').read_text())
assert result['passed']
assert json.loads((P / 'provider/build/owner.json').read_text())['passed']
assert json.loads((P / 'emission-result.json').read_text())['passed']
driver, icd = R / result['library'], R / result['ICD']
assert pin(driver) == result['library_identity']
assert pin(icd) == result['ICD_identity']
catalogue = P / 'compiler-result.json'
modules = json.loads(catalogue.read_text())['modules']
inputs = {key(p) for p in W.iterdir() if p.is_file()}
inputs |= {key(p) for p in B.rglob('*') if p.is_file()}
inputs |= {key(p) for p in (R / 'out/native-original-input-ownership/reused/prior-proof').rglob('*') if p.is_file()}
inputs |= {key(catalogue), key(driver), key(icd),
           'out/lavapipe-acqrel-ordering/ordinary/capture.py',
           'out/native-cpu-row-ownership/cold-frame/run_owned.py',
           'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
           'out/native-original-input-ownership/provider-reuse.json'}
inputs |= {key(p) for p in (P / 'provider').iterdir() if p.is_file()}
for kind in ['Endpoint', 'DopriPhase']:
    row = modules[kind]
    module = R / row['module']
    assert pin(module)['sha256'] == row['sha256']
    inputs |= {key(module), key(module.with_suffix('.spvasm'))}
tools = [Path(p) for p in ['/usr/bin/gdb', '/usr/bin/python3', '/usr/bin/nm',
                          '/usr/bin/objdump', '/bin/bash']]
inputs |= {'/etc/ld.so.cache', '/usr/bin/ldd',
           str(Path('/usr/bin/ldd').resolve(strict=True))}
for p in [*tools, driver, B / 'short_input_jit_probe', B / 'optimizer_observer.so']:
    inputs |= {str(p), str(p.resolve(strict=True))}
    text = subprocess.check_output(['/usr/bin/ldd', str(p)], text=True)
    for item in re.findall(r'(/[^\s()]+)', text):
        q = Path(item)
        if q.is_file():
            inputs |= {str(q), str(q.resolve(strict=True))}
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
actions = {}
for label, kind in [('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
    actions[label] = {
        'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
        'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'argv': ['/usr/bin/python3', '-B', str(W / 'capture.py'), kind, '{folder}'],
        'frozen_inputs': sorted(inputs), 'software_provider': True,
        'runtime_icd': str(icd), 'runtime_provider': str(driver),
        'scope': 'One original4B guard/canary ordinary called-code capture of actual new native CPU shared-op module; unchanged IEEE32/intergroup/bulk-clear and Default engine+greedy. Bind emitted input, calls, features, clear and FP restoration before considering numerical trial. No numerical/performance/frame/release acceptance.'}
(W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'actions': len(actions), 'inputs': len(inputs)}))
