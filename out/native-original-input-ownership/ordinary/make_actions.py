from pathlib import Path
import json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/native-retained-ownership'
inputs={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}|{'/usr/bin/gdb','/usr/bin/python3','/etc/ld.so.cache','out/native-original-input-ownership/compiler-result.json','out/native-original-input-ownership/provider-reuse.json'}
inputs|={str(p.relative_to(R)) for p in B.rglob('*') if p.is_file()}
inputs|={str(p.relative_to(R)) for p in (W.parent/'reused/prior-proof').rglob('*') if p.is_file()}
for p in [Path('/usr/bin/gdb'),B/'providers/candidate/libvulkan_lvp.so']:
 s=subprocess.check_output(['ldd',str(p)],text=True);inputs|={str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)',s) if Path(x).is_file()}
inputs |= {'bin/linux-gcc/native-original-input-ownership/modules/owned-endpoint.spv', 'bin/linux-gcc/native-original-input-ownership/modules/owned-transport.spv'}
a={};revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
for label,kind in [('transport','Transport'),('endpoint','Endpoint')]:
 a[label]={'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'capture.py'),kind,'{folder}'],'frozen_inputs':sorted(inputs),'software_provider':True,'runtime_icd':str(B/'providers/candidate/lvp_icd.json'),'runtime_provider':str(B/'providers/candidate/libvulkan_lvp.so'),'scope':'Actual complete '+kind+' ownership module, ordinary optimized LLVM and called JIT code. Same reused actual IEEE32 provider and public-factory guard ELF; no new arithmetic algorithm or workgroup. Guard dispatch only, no retained/scientific/benefit/frame qualification.'}
assert not (W/'actions.json').exists();(W/'actions.json').write_text(json.dumps(a,indent=2)+'\n')
