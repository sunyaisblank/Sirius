from pathlib import Path
import json,re,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;P=R/'out/lavapipe-vector-producer'
for n in ['run_owned.py','owned_guard.py']:shutil.copyfile(P/n,W/n)
p=W/'run_owned.py';s=p.read_text().replace("'VK_DRIVER_FILES':'/usr/share/vulkan/icd.d/lvp_icd.json'","'VK_DRIVER_FILES':action['runtime_icd']")
s=s.replace("provider = '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so'","provider = action['runtime_provider']")
p.write_text(s)
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
for rel in ['bin/linux-gcc/lavapipe-vector-validation','bin/linux-gcc/lavapipe-vector-staged','out/lavapipe-vector-validation-build','out/lavapipe-vector-inspect']:
 inputs|={str(p.relative_to(R)) for p in (R/rel).iterdir() if p.is_file()}
for arm in ['baseline','candidate']:
 inputs|={f'bin/linux-gcc/lavapipe-vector-staged/{arm}/libvulkan_lvp.so'}
inputs|={'out/lavapipe-vector-capture/ordinary-code-review.json','out/lavapipe-vector-link/matched-links.json','bin/linux-gcc/tests/backend/portable_binary32_reference.bin','bin/linux-gcc/tests/backend/resource-index.json','/etc/ld.so.cache','/usr/bin/python3'}
# The resource index name varies by role; use the actual existing executable-adjacent resources.
inputs={n for n in inputs if (R/n).is_file()}
for p in [R/'bin/linux-gcc/lavapipe-vector-validation/audited_backend_tests',R/'bin/linux-gcc/lavapipe-vector-validation/gpu_probe',R/'bin/linux-gcc/lavapipe-vector-staged/candidate/libvulkan_lvp.so']:
 s=subprocess.check_output(['/usr/bin/ldd',str(p)],text=True)
 inputs|={str(Path(a).resolve(strict=True)) for a in re.findall(r'(/[^\s()]+)',s) if Path(a).is_file()}
inputs|={str(p.relative_to(R)) for p in (R/'bin/linux-gcc/src/sirius/backend/retained').glob('*.spv')}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();actions={}
for case,arm in [('candidate-corpus64','candidate'),('baseline-barrier16','baseline'),('candidate-barrier16','candidate'),('control-8','candidate'),('a1','baseline'),('b1','candidate'),('b2','candidate'),('a2','baseline')]:
 corpus=('corpus' in case or 'barrier16' in case)
 prior=[] if corpus else ['candidate-corpus64','baseline-barrier16','candidate-barrier16']
 if case in ['a1','b1','b2','a2']:prior+=['control-8']+['a1','b1','b2','a2'][:['a1','b1','b2','a2'].index(case)]
 actions[case]={'revision':revision,'kind':'diagnostic' if corpus else 'test','expected_exit':0,
  'limits':{'seconds':180 if corpus else 300,'rss_kib':4194304,'reserve_kib':2097152},
  'argv':['/usr/bin/python3','-B',str(W/'run_action.py'),case,'{folder}'],'frozen_inputs':sorted(inputs),
  'software_provider':True,'runtime_icd':str(R/f'out/lavapipe-vector-inspect/{arm}-icd.json'),
  'runtime_provider':str(R/f'bin/linux-gcc/lavapipe-vector-staged/{arm}/libvulkan_lvp.so'),
  'prerequisites':prior,'scope':('Original174268 arithmetic records plus48 canary words under matched '+arm+' provider; '+('original64 fallback' if 'corpus' in case else 'cross-lane shared/atomic16')+' diagnostic only.' if corpus else 'Original C8 with exact original assertions and existing extra module/API observations.300s/4GiB/2GiB reserve, cache disabled, no debugger/preload, three original internal repeats; '+case+'. No frame/release qualification.')}
(W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps({'prepared':True,'frozen_inputs':len(inputs),'fixed_ABBA':['a1','b1','b2','a2'],'no_execution_by_preparation':True}))
