from pathlib import Path
import json,re,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64-matched'
for n in ['run_owned.py','owned_guard.py']:shutil.copyfile(R/'out/portable-checked-fp64'/n,W/n)
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
for rel in ['out/portable-checked-fp64-matched','bin/linux-gcc/portable-checked-fp64-matched']:
 inputs|={str(p.relative_to(R)) for p in (R/rel).rglob('*') if p.is_file()}
inputs|={'out/portable-checked-fp64-jit/ordinary-code-review.json','out/portable-checked-fp64/run/complete-output-agreement.json','/usr/share/vulkan/icd.d/lvp_icd.json','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/etc/ld.so.cache','/usr/bin/python3'}
inputs|={str(p.relative_to(R)) for p in (R/'bin/linux-gcc/src/sirius/backend/retained').glob('*.spv')}
inputs|={'bin/linux-gcc/portable-checked-fp64-modules/checked-'+s+'.spv' for s in ['endpoint','transport','dopri_phase']}
for exe in [B/'candidate-tests',Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')]:
 s=subprocess.check_output(['/usr/bin/ldd',str(exe)],text=True);inputs|={str(Path(a).resolve(strict=True)) for a in re.findall(r'(/[^\s()]+)',s) if Path(a).is_file()}
actions={};revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
for label,arm in [('control-8','candidate'),('a1','baseline'),('b1','candidate'),('b2','candidate'),('a2','baseline')]:
 prior=[] if label=='control-8' else ['control-8']+['a1','b1','b2','a2'][:['a1','b1','b2','a2'].index(label)]
 actions[label]={'revision':revision,'kind':'test','expected_exit':0,'arm':arm,'prerequisites':prior,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'run_case.py'),label,'{folder}'],'frozen_inputs':sorted(inputs),'software_provider':True,'scope':'Original C8 or sole fixed ABBA slot '+label+' under actual stock provider and enabled Float64. Three original internal repeats, unchanged full assertions/work/count/governor/cache and300s/4GiB/2GiB-reserve owner. Only portable normal Endpoint/Transport/DopriPhase select checked estimates; no scientific/frame/release qualification.'}
(W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps({'frozen_inputs':len(inputs),'C8_before_sole_ABBA':True,'fixed_order':['a1','b1','b2','a2']}))
