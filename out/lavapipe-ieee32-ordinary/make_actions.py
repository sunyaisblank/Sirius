from pathlib import Path
import json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-ieee32-staged'
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{'bin/linux-gcc/lavapipe-ieee32-controls/short_input_jit_probe','bin/linux-gcc/src/sirius/backend/retained/retained_transport.spv','bin/linux-gcc/lavapipe-native-fp32/float-layout.json','out/lavapipe-ieee32-validation/candidate-default/numerical-and-boundary-result.json','out/lavapipe-ieee32-validation/candidate-poison/numerical-and-boundary-result.json','/usr/bin/gdb','/usr/bin/python3','/etc/ld.so.cache'}
inputs|={str(p.relative_to(R)) for p in B.rglob('*') if p.is_file()}
for p in [Path('/usr/bin/gdb'),B/'candidate/libvulkan_lvp.so']:
    s=subprocess.check_output(['/usr/bin/ldd',str(p)],text=True);inputs|={str(Path(x).resolve(strict=True)) for x in re.findall(r'(/[^\s()]+)',s) if Path(x).is_file()}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
(W/'actions.json').write_text(json.dumps({'capture':{'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'capture.py'),'{folder}'],'frozen_inputs':sorted(inputs),'software_provider':True,'runtime_icd':str(B/'candidate/lvp_icd.json'),'runtime_provider':str(B/'candidate/libvulkan_lvp.so'),'scope':'One ordinary complete shipped native Transport code/feature binding under the actual scoped IEEE32 provider. Original passive observer and public factory ELF reused exactly. Short-input dispatch only; independent scalar corpus is separate. No performance/retained/frame acceptance.'}},indent=2)+'\n')
