from pathlib import Path
import json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-ieee32-staged'
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{str(p.relative_to(R)) for p in B.rglob('*') if p.is_file()}
inputs|={'bin/linux-gcc/lavapipe-native-fp32/gpu_probe','bin/linux-gcc/lavapipe-fpstate-contract/gpu_probe','bin/linux-gcc/lavapipe-native-fp32/native.spv','bin/linux-gcc/lavapipe-native-fp32/float-layout.json','out/lavapipe-native-fp32/input.bin','out/lavapipe-native-fp32/expected.bin','out/lavapipe-ieee32/build/owner.json','out/lavapipe-ieee32/candidate.patch','out/lavapipe-ieee32/source-admission.json','/usr/bin/gdb','/usr/bin/python3','/etc/ld.so.cache'}
for p in [Path('/usr/bin/gdb'),B/'baseline/libvulkan_lvp.so',B/'candidate/libvulkan_lvp.so']:
    s=subprocess.check_output(['/usr/bin/ldd',str(p)],text=True);inputs|={str(Path(x).resolve(strict=True)) for x in re.findall(r'(/[^\s()]+)',s) if Path(x).is_file()}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();a={}
for name,arm,poison in [('baseline-poison-v21','baseline',True),('candidate-poison','candidate',True),('candidate-default','candidate',False)]:
    a[name]={'revision':revision,'arm':arm,'poison':poison,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'run_case.py'),name,'{folder}'],'frozen_inputs':sorted(inputs),'software_provider':True,'runtime_icd':str(B/arm/'lvp_icd.json'),'runtime_provider':str(B/arm/'libvulkan_lvp.so'),'scope':('Reproduce the real normal-domain RTE32 counterexample on the exact reused unmodified LLVM21 baseline. Diagnostic completion is not numerical acceptance.' if arm=='baseline' else 'Original174268-record native arithmetic corpus under actual new supported DenormPreserve32/RTE32/SignedZeroInfNanPreserve32, with passive real worker/callback-state restoration capture. All522804 expected words exact; no retained/benefit/frame/scientific/release acceptance.')}
(W/'actions.json').write_text(json.dumps(a,indent=2)+'\n')
