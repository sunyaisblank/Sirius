"""Bind one complete-code capture per exact reused shader to the staged provider."""
from pathlib import Path
import hashlib, json, re, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W.parent
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
assert doc(P / 'build/owner.json')['passed']
result = doc(P / 'provider-build-result.json'); assert result['passed']
driver = R / result['library']; icd = R / result['ICD']
assert pin(driver) == result['library_identity'] and pin(icd) == result['ICD_identity']
compiler_path = R / 'out/native-cpu-row-ownership/compiler-result.json'
compiler = doc(compiler_path)
prior = R / 'out/native-original-input-ownership'
B = R / 'bin/linux-gcc/native-retained-ownership'
inputs = {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in B.rglob('*') if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in (prior / 'reused/prior-proof').rglob('*') if p.is_file()}
for record in compiler['modules'].values():
    module = R / record['module']; assert pin(module)['sha256'] == record['sha256']
    inputs |= {str(module.relative_to(R)), str(module.with_suffix('.spvasm').relative_to(R))}
inputs |= {str(p.relative_to(R)) for p in [driver, icd, compiler_path,
    P / 'build/owner.json', P / 'build/inputs-before.json', P / 'build/inputs-after.json',
    P / 'provider-build-result.json', P / 'review-admission.json', P / 'source-contract.json',
    P / 'reused-row-products.json', R / 'out/lavapipe-acqrel-ordering/ordinary/capture.py',
    prior / 'provider-reuse.json']}
inputs |= {'/usr/bin/gdb', '/usr/bin/python3', '/usr/bin/objdump', '/usr/bin/nm', '/etc/ld.so.cache',
           'out/native-cpu-row-ownership/cold-frame/run_owned.py'}
for obj in [Path('/usr/bin/gdb'), Path('/usr/bin/nm'), driver, B / 'short_input_jit_probe']:
    text = subprocess.check_output(['ldd', str(obj)], text=True)
    inputs |= {str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)', text) if Path(x).is_file()}
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
actions = {}
for label, kind in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
    actions[label] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
        'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'argv': ['/usr/bin/python3', '-B', str(W / 'capture.py'), kind, '{folder}'],
        'frozen_inputs': sorted(inputs), 'software_provider': True,
        'runtime_icd': str(icd), 'runtime_provider': str(driver),
        'scope': 'One actual ordinary complete '+kind+' code capture from the staged intergroup SIMD provider and exact preserved bfa095ba row-owner module. Original four-byte guard/canary dispatch, cache policy, factory features and callback FP-state observation. No C++/shader/probe compilation, numerical/benefit/frame/science/release acceptance. Actual divergent group IDs, original grid/block arguments and pre-prologue tail masking must be reviewed before C8.'}
target = W / 'actions.json'
if target.exists():
    # A completed diagnostic is never repeated to add an unused tool seal.
    # Preserve its exact original action document and source revision.
    assert doc(W / 'transport/owner.json')['passed']
    assert not (W / 'endpoint').exists() and not (W / 'dopri_phase').exists()
    original = W / 'transport/original-actions.json'; assert not original.exists()
    original.write_bytes(target.read_bytes())
    actions['transport'] = doc(target)['transport']
target.write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'new_provider': pin(driver), 'complete_module_captures_planned': len(actions),
                  'new_CPP_probe_shader_compilations': 0}))
