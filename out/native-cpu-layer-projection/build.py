"""Compile one exact DP layer projection and inspect its local SSA boundary."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,shlex,subprocess,sys

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking'
B=R/'bin/linux-gcc/native-cpu-strict-ir-linking';D=B/'toolchain/bin'
I=R/'bin/linux-gcc/native-cpu-local-effects/dp.cleanup.bc';O=R/'bin/linux-gcc/native-cpu-layer-projection'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def doc(p):return json.loads(p.read_text())


def write(p,v):
 with p.open('x') as s:json.dump(v,s,indent=2);s.write('\n')


phase,=sys.argv[1:]
if phase=='prepare':
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 assert not O.exists();O.mkdir()
 config=shlex.split(subprocess.check_output([str(D/'llvm-config'),'--cxxflags'],text=True))
 pipeline='always-inline,function(sroa,mem2reg,instcombine,simplifycfg,early-cse<memssa>,dse),globaldce,verify'
 commands=dict(plugin=[str(D/'clang++'),*config,'-O2','-fPIC','-shared',str(V/'layer_projection.cpp'),str(B/'sha256-only.a'),'-Wl,--exclude-libs,ALL','-o',str(O/'layer_projection.so')],
  lower=[str(D/'opt'),'-load-pass-plugin='+str(O/'layer_projection.so'),'-passes=sirius-project-dp,verify','-verify-each',str(I),'-o',str(O/'dp.bc')],
  cleanup=[str(D/'opt'),'-passes='+pipeline,'-verify-each','-mcpu=skylake-avx512',str(O/'dp.bc'),'-o',str(O/'dp.cleanup.bc')],
  native=[str(D/'llc'),'-O2','-filetype=obj','-relocation-model=pic','-mcpu=skylake-avx512',str(O/'dp.cleanup.bc'),'-o',str(O/'dp.o')])
 outputs=dict(plugin='layer_projection.so',lower='dp.bc',cleanup='dp.cleanup.bc',native='dp.o')
 paths=[*W.glob('*.py'),V/'layer_projection.cpp',V/'layer_projection_contract.json',I,B/'sha256-only.a',V/'toolchain.json',
  R/'out/native-cpu-local-effects/dp-cleanup/owner.json',R/'out/native-cpu-local-effects/dp-cleanup-result.json']
 for name in ['clang++','opt','llvm-config','llc','llvm-dis']:
  tool=doc(V/'toolchain.json')['tools'][name];p=R/tool['resolved'];assert pin(p)=={k:tool[k] for k in ['bytes','sha256']};paths.append(p)
 inputs={str(p.relative_to(R)):pin(p) for p in paths}
 write(W/'plan.json',dict(commands=commands,outputs=outputs,inputs=inputs,cleanup_pipeline=pipeline,scope='First DP-only local-slot projection gate; no numerical/workload/adoption qualification.'))
 actions={};prior=dict(plugin=None,lower='plugin',cleanup='lower',native='cleanup')
 for name in commands:
  extra=[] if not prior[name] else [str((O/outputs[prior[name]]).relative_to(R))]
  actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
   limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=[sys.executable,'-B',str(Path(__file__).resolve()),name],
   frozen_inputs=sorted(inputs)+[str((W/'plan.json').relative_to(R))]+extra,
   scope='One DP projection compiler gate: initialized private per-slot locals, exact strict FP inventory/math calls/old raw bounds, fixed compact interface and original work limits. No frame, benchmark, numerical/provider/scientific/hardware/release acceptance.')
 write(W/'actions.json',actions);print(json.dumps(dict(prepared=list(commands),runtime_run=False)))
else:
 plan=doc(W/'plan.json');assert phase in plan['commands']
 for name,expected in plan['inputs'].items():assert pin(R/name)==expected,name
 prior=dict(plugin=None,lower='plugin',cleanup='lower',native='cleanup')[phase]
 if prior:
  owner=doc(W/prior/'owner.json');assert owner['passed'] and owner['status']=='terminal'
  product=doc(W/(prior+'-result.json'))['output'];assert pin(R/product['path'])=={k:product[k] for k in ['bytes','sha256']}
 output=O/plan['outputs'][phase];assert not output.exists()
 subprocess.run(plan['commands'][phase],check=True)
 result=dict(phase=phase,output=dict(path=str(output.relative_to(R)),**pin(output)),argv=plan['commands'][phase])
 if output.suffix=='.bc':
  ll=output.with_suffix('.ll');subprocess.run([str(D/'llvm-dis'),str(output),'-o',str(ll)],check=True)
  text=ll.read_text();definitions=re.findall(r'^define .*$',text,re.M);attrs=dict(re.findall(r'^attributes #(\d+) = (.*)$',text,re.M))
  assert all('strictfp' in attrs[re.search(r'#(\d+)',line)[1]] for line in definitions)
  assert not re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',text,re.M)
  cores=list(re.finditer(r'^define [^\n]*@(sirius_projected_layer_\d+)\([^\n]*\{\n(.*?)^}',text,re.M|re.S));assert len(cores)==41
  core_calls=Counter(re.findall(r'\bcall\s+.*?@(sirius_local_value_[^ (]+)\(',text));assert sum(core_calls.values())==1581
  inventory=Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(',text))
  result.update(ir=dict(path=str(ll.relative_to(R)),**pin(ll)),definitions=len(definitions),core_calls=dict(core_calls),constrained_calls=dict(inventory),
   projected_cores=[dict(name=m[1],local_allocas=len(re.findall(r'^\s+%project.local.[^=]+ = alloca',m[2],re.M)),
    masked_calls=len(re.findall(r'@llvm.masked.(?:load|store)',m[2])),ordinary_loads=len(re.findall(r' = load ',m[2])),ordinary_stores=len(re.findall(r'^\s+store ',m[2],re.M))) for m in cores])
  if phase=='cleanup':
   assert all(r['local_allocas']==0 and r['masked_calls']==0 for r in result['projected_cores']),result['projected_cores']
   assert inventory==Counter(doc(W/'lower-result.json')['constrained_calls'])
 elif phase=='native':
  asm=subprocess.check_output(['objdump','-drw',str(output)],text=True)
  instructions=re.findall(r'^\s*[0-9a-f]+:\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9.]*)',asm,re.M)
  assert instructions and not [x for x in instructions if x.startswith(('vfm','vfnm'))]
  refs=Counter(re.findall(r'R_X86_64_PLT32\s+(sirius_local_value_\w+)-0x4',asm));assert sum(refs.values())==1581
  result.update(instructions=len(instructions),core_relocations=dict(refs),emitted_fused_instructions=[])
 write(W/(phase+'-result.json'),result);print(json.dumps({k:result[k] for k in ['phase','output','definitions','instructions'] if k in result}))
