"""Prepare isolated original stock-route controls; no compiler/device call."""
from pathlib import Path
import hashlib,json,shlex,shutil,subprocess
r=Path.cwd();w=Path(__file__).resolve().parent;b=r/'bin/linux-gcc/portable-integer-reciprocal';phase=w/'original-controls';phase.mkdir()
old=r/'attestations/software-vulkan/dcd1368/sealed-stage-seven-owner/files/out/sealed-stage-seven-owner';saved=old/'original-controls';rev=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
for mode in ('baseline','candidate'):
 target=phase/mode;target.mkdir();shutil.copyfile(saved/'baseline/retained_compute_test.cpp',target/'retained_compute_test.cpp')
 (target/'retained_dopri_test.cpp').write_bytes(subprocess.check_output(['git','show','9fbd0e8:tests/backend/retained_dopri_test.cpp']))
assert (phase/'baseline/retained_compute_test.cpp').read_bytes()==(phase/'candidate/retained_compute_test.cpp').read_bytes()
inc=phase/'baseline/include/sirius/backend';inc.mkdir(parents=True);(inc/'retained_compute.h').write_bytes(subprocess.check_output(['git','show','9fbd0e8:src/sirius/backend/retained_compute.h']))
for name in ('run_owned.py','owned_guard.py'):shutil.copyfile(w/name,phase/name)
cc=json.loads((r/'bin/linux-gcc/compile_commands.json').read_text());actions={};lim={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
compiler=['/usr/bin/g++-14','/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus','/usr/libexec/gcc/x86_64-linux-gnu/14/collect2','/usr/bin/x86_64-linux-gnu-ld.bfd',r/'bin/linux-gcc/compile_commands.json']
def add(name,argv,inputs,kind='diagnostic',soft=False,scope='Strict original-control private preparation; no canonical configure/build'):
 actions[name]={'revision':rev,'kind':kind,'expected_exit':0,'argv':list(map(str,argv)),'limits':lim,'frozen_inputs':list(map(str,[*compiler,*inputs])),'scope':scope}
 if soft:actions[name]['software_provider']=True
for mode in ('baseline','candidate'):
 for file,tag in [('retained_compute_test.cpp','observer'),('retained_dopri_test.cpp','dopri')]:
  command=shlex.split(next(x['command'] for x in cc if x['file'].endswith('/'+file)))
  obj=b/(mode+'-'+tag+'.o');command[command.index('-o')+1]=str(obj);command[command.index('-c')+1]=str(phase/mode/file)
  command.insert(1,'-I'+str(phase/'baseline/include' if mode=='baseline' else b/'products'))
  command+=['-MMD','-MF',str(obj)+'.d']
  add('compile-'+mode+'-'+tag,command,[phase/mode/file,*((phase/'baseline/include').rglob('*.h') if mode=='baseline' else [b/'products/retained_kernels.h'])],kind='compiler')
for file,tag in [('retained_compute.cpp','compute'),('retained_integrator.cpp','integrator'),('retained_trace_executor.cpp','trace')]:
 command=shlex.split(next(x['command'] for x in cc if x['file'].endswith('/'+file)));obj=b/(tag+'.o');command[command.index('-o')+1]=str(obj);command.insert(1,'-I'+str(b/'products'));command+=['-MMD','-MF',str(obj)+'.d'];add('compile-'+tag,command,[b/'products/retained_kernels.h'],kind='compiler')
libs=['src/sirius/backend/libsirius_backend_cpu.a','lib/libgtest_main.a','src/sirius/backend/libsirius_backend.a','src/sirius/render/libsirius_render.a','lib/libgtest.a','src/sirius/backend/libsirius_backend_cpu.a','src/sirius/backend/libsirius_backend.a','/usr/lib/x86_64-linux-gnu/libvulkan.so','-ldl','src/sirius/core/libsirius_core.a','src/sirius/libsirius_base.a']
libs=[str(r/'bin/linux-gcc'/p) if not p.startswith(('/','-')) else p for p in libs]
for mode in ('baseline','candidate'):
 objects=[b/(mode+'-observer.o'),b/(mode+'-dopri.o')]
 if mode=='candidate':objects += [b/(p+'.o') for p in ('compute','integrator','trace')]
 objects += [b/'vulkan.o'];exe=b/('original_'+mode+'_controls')
 add('link-'+mode,['/usr/bin/g++-14','-O3','-DNDEBUG',*objects,*libs,'-pthread','-Wl,-Map='+str(b/(mode+'.map')),'-o',exe],[*objects,*[p for p in libs if not p.startswith('-')]],scope='Private stock-route link. Baseline9fbd source/products and candidatecurrent products share exact current Vulkan object/provider; candidate overrides all old Compute/integrator/trace archive members.')
oldactions=json.loads((saved/'actions.json').read_text())
for n in range(1,9):
 v=oldactions['control-'+str(n)];v['revision']=rev;v['argv'][0]=str(b/'original_candidate_controls');v['frozen_inputs']=[str(b/'original_candidate_controls'),'/usr/lib/x86_64-linux-gnu/libvulkan.so','/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/share/vulkan/icd.d/lvp_icd.json'];v['scope']='Unchanged original finite scientific/source/ordered-stream/refusal/recovery filter on changed shared portable arithmetic; not full science/release qualification';actions['control-'+str(n)]=v
(phase/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(phase/'source-bindings.json').write_text(json.dumps({'revision':rev,'observer_origin':str(saved/'baseline/retained_compute_test.cpp'),'observer_sha256':sha(saved/'baseline/retained_compute_test.cpp'),'baseline_and_candidate_observer_identical':True,'dopri_source_revision':'9fbd0e8','default_stock_route':True,'owner_opt_in':False,'all_compilation_limits':lim},indent=2)+'\n')
print(json.dumps({'actions':list(actions),'scope':'preparation only; no compiler or device dispatch'}))
