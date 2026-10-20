from pathlib import Path
import json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-fpstate-contract'
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
for row in json.loads((W/'commands.json').read_text()):
    for arg in row['argv']:
        p=Path(arg);p=p if p.is_absolute() else Path(row['cwd'])/p
        if p.is_file():inputs.add(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()))
inputs|={'bin/linux-gcc/lavapipe-fpstate-contract/'+n for n in ['native.spv','native.spvasm','float-layout.json']}
inputs|={'/usr/bin/g++-14','/usr/bin/gdb','/usr/bin/python3','/usr/share/vulkan/icd.d/lvp_icd.json','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275','/etc/ld.so.cache','attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so'}
for p in [Path('/usr/bin/gdb'),Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')]:
    s=subprocess.check_output(['/usr/bin/ldd',str(p)],text=True);inputs|={str(Path(x).resolve(strict=True)) for x in re.findall(r'(/[^\s()]+)',s) if Path(x).is_file()}
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();a={}
for name in ['build','capture']:
    f=inputs|({'bin/linux-gcc/lavapipe-fpstate-contract/gpu_probe'} if name=='capture' else set())
    a[name]={'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/('build.py' if name=='build' else 'run_capture.py'))]+([] if name=='build' else ['{folder}']),'frozen_inputs':sorted(f),'software_provider':name=='capture','scope':('Compile only the original strict factory probe with scoped hostile host rounding/denormal state; exact original native SPV and oracle reused.' if name=='build' else 'One changed-input contract observation: caller uses round-up/FTZ/DAZ while actual supported shader RTE32 is explicitly requested. Passive worker MXCSR/property and original-output capture; no driver mode/property override. Owner PASS means diagnostic completed; mismatches remain numerical failures. No rendering/benefit/feature qualification.')}
(W/'actions.json').write_text(json.dumps(a,indent=2)+'\n')
