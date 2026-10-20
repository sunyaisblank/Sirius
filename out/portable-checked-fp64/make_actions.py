from pathlib import Path
import json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64'
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{'bin/linux-gcc/tests/backend/portable_binary32_reference.bin'}
for row in json.loads((W/'commands.json').read_text()):
    for arg in row['argv']:
        p=Path(arg);p=p if p.is_absolute() else Path(row['cwd'])/p
        if p.is_file():inputs.add(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()))
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();a={}
for mode in ['build','host','run']:
    f=inputs|({str((B/name).relative_to(R)) for name in ['check_host','gpu_probe.o','gpu_probe','corpus.spv','corpus.spvasm']} if mode!='build' else set())
    if mode=='run':f|={'/usr/share/vulkan/icd.d/lvp_icd.json','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1','/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275'}
    a[mode]={'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':120,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'step.py'),mode,'{folder}'],'frozen_inputs':sorted(f),'software_provider':mode=='run','scope':{'build':'Compile private exact-verified FP64 estimates with strict existing C++ probe and original scalar shader; validate actual capability/types and preserve all integer fallbacks.','host':'One host check:720288 deliberately good/bad seeds against independent integer quotients and binary-search roots, plus original174268-record corpus. UBSan and float-cast-overflow enabled; no Vulkan.','run':'One original174268-record scalar dispatch with actual supported/enabled Float64 public factory. All522804 original expected words complemented before execution and compared exactly. No integrated, benefit, frame or release admission.'}[mode]}
(W/'actions.json').write_text(json.dumps(a,indent=2)+'\n')
