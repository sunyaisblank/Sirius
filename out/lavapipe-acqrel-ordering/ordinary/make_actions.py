"""Freeze the real new provider with exact previously admitted complete modules."""
from pathlib import Path
import hashlib, json, re, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W.parent
B = R / 'bin/linux-gcc/native-retained-ownership'
driver = R / 'bin/linux-gcc/lavapipe-acqrel-ordering/src/gallium/targets/lavapipe/libvulkan_lvp.so'
assert json.loads((P / 'build/owner.json').read_text())['passed'] and driver.is_file()
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
icd = P / 'provider-icd.json'; assert not icd.exists()
icd.write_text(json.dumps({'file_format_version': '1.0.0',
    'ICD': {'library_path': str(driver), 'api_version': '1.4.305'}}, indent=2) + '\n')
prior = R / 'out/native-original-input-ownership'
compiler = json.loads((prior / 'compiler-result.json').read_text())
modules = [R / x['module'] for x in compiler['modules'].values()]
for kind, record in compiler['modules'].items():
    assert pin(R / record['module'])['sha256'] == record['sha256']
inputs = {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in B.rglob('*') if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in (prior / 'reused/prior-proof').rglob('*') if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in modules}
inputs |= {str(p.relative_to(R)) for p in [driver, icd, P / 'build/owner.json',
    P / 'build/inputs-before.json', P / 'build/inputs-after.json', P / 'source-admission.json',
    P / 'source-contract.json', P / 'configuration-command.json', prior / 'compiler-result.json',
    prior / 'provider-reuse.json']}
inputs |= {'/usr/bin/gdb', '/usr/bin/python3', '/usr/bin/objdump', '/etc/ld.so.cache'}
for obj in [Path('/usr/bin/gdb'), driver, B / 'short_input_jit_probe']:
    text = subprocess.check_output(['ldd', str(obj)], text=True)
    inputs |= {str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)', text) if Path(x).is_file()}
actions = {}; revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
for label, kind in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
    actions[label] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
        'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'argv': ['/usr/bin/python3', '-B', str(W / 'capture.py'), kind, '{folder}'],
        'frozen_inputs': sorted(inputs), 'software_provider': True,
        'runtime_icd': str(icd), 'runtime_provider': str(driver),
        'scope': 'Actual exact previous complete ' + kind + ' module with new private acquire-release provider; ordinary optimized LLVM/called code, actual factory precision and callback restoration. Four-byte guard only. Preserve cooperative64-lane workgroup, rendezvous, atomics, programs and arithmetic. No additional shader/C++ compilation, retained numerical/benefit/frame/science/release acceptance.'}
assert not (W / 'actions.json').exists()
(W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
(P / 'provider-identity.json').write_text(json.dumps({'library': str(driver.relative_to(R)),
    'identity': pin(driver), 'ICD': str(icd.relative_to(R)), 'ICD_identity': pin(icd),
    'unchanged_shaders': {k: {'module': v['module'], 'sha256': v['sha256']} for k, v in compiler['modules'].items()},
    'original_precision_proof': '6a6bb1ca48693e290999a0aab8f305cc26d91262; unchanged exact5-file IEEE32 prerequisite, composed with new actual provider observations.',
    'build_owner': pin(P / 'build/owner.json'), 'provider_installed': False,
    'provider_or_source_adopted': False}, indent=2) + '\n')
print(json.dumps({'new_provider': pin(driver), 'unchanged_complete_modules': len(modules),
                  'probe_compilations': 0, 'shader_compilations': 0}))
