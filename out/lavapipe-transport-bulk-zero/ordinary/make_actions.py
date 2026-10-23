"""Materialize one finite code screen, with exact actual probe/provider and passive authorities."""
from pathlib import Path
import hashlib,json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;P=W.parent
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def key(p):return str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
assert not (W/'actions.json').exists() and not (W/'transport').exists()
assert not any(p.is_symlink() and p.is_dir() for p in W.iterdir())
result=doc(P/'provider-build-result.json');assert result['passed'] and doc(P/'build/owner.json')['passed']
driver=R/result['library'];icd=R/result['ICD']
assert pin(driver)==result['library_identity'] and pin(icd)==result['ICD_identity']
catalogue=R/'out/native-cpu-row-ownership/compiler-result.json';c=doc(catalogue)
module=R/c['modules']['Transport']['module'];assert pin(module)=={k:c['modules']['Transport'][k] for k in ['bytes','sha256']}
assert hashlib.sha1(module.read_bytes()).hexdigest()==doc(P/'source-contract.json')['Transport']['SPIRV_SHA1']
B=R/'bin/linux-gcc/native-retained-ownership'
inputs={key(p) for p in W.iterdir() if p.is_file()}
inputs|={key(p) for p in B.rglob('*') if p.is_file()}
prior=R/'out/native-original-input-ownership'
inputs|={key(p) for p in (prior/'reused/prior-proof').rglob('*') if p.is_file()}
proofs=[P/'provider-build-result.json',P/'source-contract.json',P/'bulk-zero.patch',P/'review-admission.json',
 P/'build/owner.json',P/'build/inputs-before.json',P/'build/inputs-after.json',catalogue,module,module.with_suffix('.spvasm'),
 prior/'provider-reuse.json',R/'out/lavapipe-acqrel-ordering/ordinary/capture.py',
 R/'out/lavapipe-owner-bulk-zero/ordinary/capture.py',R/'out/lavapipe-owner-bulk-zero/ordinary/symbol_bindings.py',
 R/'out/lavapipe-owner-bulk-zero/ordinary/read_code.py',R/'out/native-cpu-row-ownership/cold-frame/run_owned.py',
 R/'out/lavapipe-intergroup-row-owner/ordinary/ordinary-code-review.json',driver,icd]
root=R/'out/lavapipe-intergroup-row-owner/ordinary/transport'
v=doc(root/'capture-verification.json');assert v['passed'] and v['actual_complete_module']['sha256']==c['modules']['Transport']['sha256']
ir=root/v['ordinary_ir_export']['path'];assert pin(ir)=={k:v['ordinary_ir_export'][k] for k in ['bytes','sha256']}
proofs += [root/'capture-verification.json',root/'bindings.json',root/'owner.json',ir]
proofs += [R/r['candidate']['path'] for r in doc(P/'source-contract.json')['modified_files'].values()]
inputs|={key(p) for p in proofs}
inputs|={'/etc/ld.so.cache','/usr/bin/ldd',str(Path('/usr/bin/ldd').resolve(strict=True))}
for obj in [*[Path(x) for x in ['/usr/bin/gdb','/usr/bin/python3','/usr/bin/nm','/usr/bin/objdump','/bin/bash']],
 driver,B/'short_input_jit_probe',B/'optimizer_observer.so']:
 inputs|={str(obj),str(obj.resolve(strict=True))}
 output=subprocess.check_output(['/usr/bin/ldd',str(obj)],text=True)
 for name in re.findall(r'(/[^\s()]+)',output):
  p=Path(name)
  if p.is_file():inputs|={str(p),str(p.resolve(strict=True))}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert not subprocess.check_output(['git','status','--porcelain']).strip()
a={'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},
 'argv':['/usr/bin/python3','-B',str(W/'capture.py'),'Transport','{folder}'],'frozen_inputs':sorted(inputs),
 'software_provider':True,'runtime_icd':str(icd),'runtime_provider':str(driver),
 'scope':'One original4B guard/canary, exact full immutable row-T module and newly staged provider. Actual complete bulkclear/Default engine/remaining-IR gate only. Retained math is not executed; no numerical/performance/frame/science/release claim.'}
(W/'actions.json').write_text(json.dumps({'transport':a},indent=2)+'\n')
print(json.dumps({'passed':True,'capture_actions':1,'explicit_input_selectors':len(inputs),'provider':pin(driver)}))
