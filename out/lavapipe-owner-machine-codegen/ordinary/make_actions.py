"""ROOT-only materialization of three finite captures; never executes a probe."""
from pathlib import Path
import hashlib, json, os, re, subprocess
R = Path.cwd(); W = Path(__file__).resolve().parent; P = W.parent
def doc(p): return json.loads(p.read_text())
def pin(p):
    with p.open('rb') as stream: digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': p.stat().st_size, 'sha256': digest}
def key(p): return str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
assert not (W/'actions.json').exists()
assert not any((W/x).exists() for x in ['transport','endpoint','dopri_phase'])
assert not any(p.is_symlink() and p.is_dir() for p in W.iterdir())
result = doc(P/'provider-build-result.json'); assert result['passed']
assert doc(P/'build/owner.json')['passed'] and doc(P/'review-admission.json')['passed']
driver = R/result['library']; icd = R/result['ICD']
assert pin(driver) == result['library_identity'] and pin(icd) == result['ICD_identity']
host = R/'bin/linux-gcc/native-cooperative-transport-intergroup/matched/audited_backend_tests'
assert pin(host)['sha256'] == 'f958a60563b64c79d3f2fb6da2815bde0effdbbe49f121034cfd212ab17b27ff'
catalogue = R/'out/native-cooperative-transport-intergroup/compiler-result.json'
compiler = doc(catalogue)
B = R/'bin/linux-gcc/native-retained-ownership'
inputs = {key(p) for p in W.iterdir() if p.is_file()}
inputs |= {key(p) for p in B.rglob('*') if p.is_file()}
prior = R/'out/native-original-input-ownership'
inputs |= {key(p) for p in (prior/'reused/prior-proof').rglob('*') if p.is_file()}
proofs = [P/'provider-build-result.json', P/'source-contract.json', P/'machine-codegen.patch',
    P/'review-admission.json', P/'build/owner.json', P/'build/inputs-before.json',
    P/'build/inputs-after.json', catalogue, host, prior/'provider-reuse.json',
    R/'out/lavapipe-acqrel-ordering/ordinary/capture.py',
    R/'out/native-cpu-row-ownership/cold-frame/run_owned.py',
    R/'out/lavapipe-intergroup-row-owner/ordinary/ordinary-code-review.json',
    prior/'ordinary/ordinary-code-review.json', driver, icd]
for label, kind in [('transport','Transport'),('endpoint','Endpoint'),('dopri_phase','DopriPhase')]:
    module = R/compiler['modules'][kind]['module']
    assert pin(module)['sha256'] == compiler['modules'][kind]['sha256']
    proofs += [module, module.with_suffix('.spvasm')]
    root = prior/'ordinary' if kind == 'Transport' else R/'out/lavapipe-intergroup-row-owner/ordinary'
    v = doc(root/label/'capture-verification.json')
    assert v['passed'] and v['actual_complete_module']['sha256'] == compiler['modules'][kind]['sha256']
    ir = root/label/v['ordinary_ir_export']['path']
    assert pin(ir) == {k:v['ordinary_ir_export'][k] for k in ['bytes','sha256']}
    proofs += [root/label/'capture-verification.json', root/label/'bindings.json',
               root/label/'owner.json', ir]
for record in doc(P/'source-contract.json')['modified_files'].values():
    proofs.append(R/record['candidate']['path'])
inputs |= {key(p) for p in proofs}
# These are actual executable paths, not version observations. Freeze resolved
# aliases and fresh actual ldd closures before the first capture, including nm.
tools = [Path(x) for x in ['/usr/bin/gdb','/usr/bin/python3','/usr/bin/nm','/usr/bin/objdump','/bin/bash']]
inputs |= {'/etc/ld.so.cache', '/usr/bin/ldd', str(Path('/usr/bin/ldd').resolve(strict=True))}
for obj in [*tools, driver, B/'short_input_jit_probe', B/'optimizer_observer.so']:
    inputs |= {str(obj), str(obj.resolve(strict=True))}
    output = subprocess.check_output(['/usr/bin/ldd', str(obj)], text=True)
    for name in re.findall(r'(/[^\s()]+)', output):
        p = Path(name)
        if p.is_file(): inputs |= {str(p), str(p.resolve(strict=True))}
revision = subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip()
assert not subprocess.check_output(['git','status','--porcelain']).strip()
actions = {}
for label, kind in [('transport','Transport'),('endpoint','Endpoint'),('dopri_phase','DopriPhase')]:
    actions[label] = {'revision':revision, 'kind':'diagnostic', 'expected_exit':0,
        'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},
        'argv':['/usr/bin/python3','-B',str(W/'capture.py'),kind,'{folder}'],
        'frozen_inputs':sorted(inputs), 'software_provider':True,
        'runtime_icd':str(icd), 'runtime_provider':str(driver),
        'scope':'One original 4B guard/canary complete-code capture; exact mixed modules and new staged provider. Engine/factory policy isolation plus E/DP optimized-IR identity only. Frozen probe is distinct from pinned intended f958 host. No retained numerical, benefit, frame or qualification conclusion.'}
(W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
print(json.dumps({'passed':True,'capture_actions':3,'explicit_input_selectors':len(inputs),'provider':pin(driver)}))
