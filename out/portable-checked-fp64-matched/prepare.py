from pathlib import Path
import json,shlex,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;C=R/'bin/linux-gcc';B=C/'portable-checked-fp64-matched'
rows=json.loads((C/'compile_commands.json').read_text());commands=[]
for suffix,source,output,extra in [('/backend/retained_compute.cpp',W/'retained_compute.cpp',B/'candidate/retained_compute.cpp.o',['-I'+str(W/'headers')])]:
 row,=[x for x in rows if x['file'].endswith(suffix)];args=shlex.split(row['command']);args[args.index('-o')+1]=str(output);args[args.index('-c')+1]=str(source);args[1:1]=extra;commands.append({'argv':args,'cwd':row['directory']})
lines=subprocess.check_output(['ninja','-C',str(C),'-t','commands','tests/backend/sirius_backend_tests'],text=True).splitlines();line,=[x for x in lines if ' -o tests/backend/sirius_backend_tests ' in x]
a=shlex.split(line);start=a.index('/usr/bin/g++-14');end=a.index('&&',start);a=a[start:end]
old='tests/backend/CMakeFiles/sirius_backend_tests.dir/retained_compute_test.cpp.o';assert a.count(old)==1
a[a.index(old)]=str(B/'audited-observer.o');links={}
for arm in ['baseline','candidate']:
 args=list(a);args[args.index('-o')+1]=str(B/(arm+'-tests'))
 if arm=='candidate':
  assert args.count('src/sirius/backend/libsirius_backend.a')==2
  args=[str(B/'candidate/libsirius_backend.a') if x=='src/sirius/backend/libsirius_backend.a' else x for x in args]
 links[arm]={'argv':args,'cwd':str(C)}
(W/'commands.json').write_text(json.dumps({'compile':commands,'link':links},indent=2)+'\n')
for n in ['run_owned.py','owned_guard.py']:shutil.copyfile(R/'out/portable-checked-fp64'/n,W/n)
inputs={str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()}
for row in commands+list(links.values()):
 for arg in row['argv']:
  p=Path(arg);p=p if p.is_absolute() else Path(row['cwd'])/p
  if p.is_file():inputs.add(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()))
inputs|={str(p.relative_to(R)) for p in (C/'src/sirius/backend/retained').iterdir() if p.is_file()}
inputs|={'bin/linux-gcc/portable-checked-fp64-modules/checked-endpoint.spv','bin/linux-gcc/portable-checked-fp64-modules/checked-transport.spv','bin/linux-gcc/portable-checked-fp64-modules/checked-dopri_phase.spv','/usr/bin/ar','out/portable-checked-fp64-jit/ordinary-code-review.json','out/portable-checked-fp64/run/complete-output-agreement.json'}
a={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'kind':'diagnostic','expected_exit':0,'limits':{'seconds':180,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'build.py')],'frozen_inputs':sorted(inputs),'software_provider':False,'scope':'Compile only one candidate three-stage selector/object; reuse the byte-exact original audited observer and all other GCC objects. Preserve all35 old arrays, add three uniquely named checked-FP64 modules, verify exactly one archive member changes. No Vulkan dispatch or integrated/benefit qualification.'}
(W/'actions.json').write_text(json.dumps({'build':a},indent=2)+'\n');print(json.dumps({'frozen_inputs':len(inputs),'new_cpp_units':1,'observer_compilations':0}))
