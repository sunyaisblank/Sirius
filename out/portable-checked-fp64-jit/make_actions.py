from pathlib import Path
import json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{'bin/linux-gcc/portable-checked-fp64-jit/short_input_jit_probe','bin/linux-gcc/portable-checked-fp64-modules/checked-endpoint.spv','out/portable-checked-fp64-modules/compiler-result.json','out/portable-checked-fp64/run/complete-output-agreement.json','/usr/bin/gdb','/usr/bin/python3','/usr/share/vulkan/icd.d/lvp_icd.json','/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/etc/ld.so.cache','attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so'}
for p in [Path('/usr/bin/gdb'),Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')]:
    output=subprocess.check_output(['/usr/bin/ldd',str(p)],text=True)
    inputs|={str(Path(x).resolve(strict=True)) for x in re.findall(r'(/[^\s()]+)',output) if Path(x).is_file()}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
(W/'actions.json').write_text(json.dumps({'capture-candidate':{'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'run_capture.py'),'{folder}','candidate'],'frozen_inputs':sorted(inputs),'software_provider':True,'scope':'One ordinary complete checked-FP64 Endpoint JIT using byte-exact preserved public factory probe and original passive observer. Actual Float64 query/create/module/code/dispatch binding required; 4B guard input is not retained arithmetic or benefit qualification.'}},indent=2)+'\n')
