"""Refine proved local ordinary effects; inspect first DP cleanup before adoption."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,shlex,subprocess,sys

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking'
I=R/'bin/linux-gcc/native-cpu-register-value';O=R/'bin/linux-gcc/native-cpu-local-effects'
D=R/'bin/linux-gcc/native-cpu-strict-ir-linking/toolchain/bin';PREFIX='sirius_local_value_'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def doc(p):return json.loads(p.read_text())


def write(p,v):
 with p.open('x') as s:json.dump(v,s,indent=2);s.write('\n')


phase,=sys.argv[1:]
if phase=='prepare':
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 assert not O.exists();O.mkdir()
 review=doc(V/'private_finite_closure.json')
 config=shlex.split(subprocess.check_output([str(D/'llvm-config'),'--cxxflags'],text=True))
 commands=dict(plugin=[str(D/'clang++'),*config,'-O2','-fPIC','-shared',str(V/'local_value_effects.cpp'),str(D.parent.parent/'sha256-only.a'),'-Wl,--exclude-libs,ALL','-o',str(O/'local_value_effects.so')])
 outputs=dict(plugin='local_value_effects.so')
 for stage in ['math','dp']:
  name=stage+'-lower';outputs[name]=stage+'.bc'
  commands[name]=[str(D/'opt'),'-load-pass-plugin='+str(O/'local_value_effects.so'),'-passes=sirius-local-'+stage+',verify','-verify-each',str(I/(stage+'.bc')),'-o',str(O/outputs[name])]
 outputs['dp-cleanup']='dp.cleanup.bc'
 pipeline='function(sroa,early-cse<memssa>,instcombine,simplifycfg,dse,early-cse<memssa>),globaldce,verify'
 commands['dp-cleanup']=[str(D/'opt'),'-passes='+pipeline,'-verify-each','-mcpu=skylake-avx512',str(O/'dp.bc'),'-o',str(O/'dp.cleanup.bc')]
 for stage in ['math','dp']:
  name=stage+'-native';outputs[name]=stage+'.o'
  commands[name]=[str(D/'llc'),'-O2','-filetype=obj','-relocation-model=pic','-mcpu='+('x86-64-v4' if stage=='math' else 'skylake-avx512'),str(O/('dp.cleanup.bc' if stage=='dp' else 'math.bc')),'-o',str(O/outputs[name])]
 paths=[*W.glob('*.py'),V/'local_value_effects.cpp',V/'private_finite_closure.json',V/'toolchain.json',D.parent.parent/'sha256-only.a']
 paths+=[I/(stage+'.bc') for stage in ['math','dp','transport','endpoint']]
 for name in ['clang++','opt','llvm-config','llc','llvm-dis']:
  tool=doc(V/'toolchain.json')['tools'][name];p=R/tool['resolved'];assert pin(p)=={k:tool[k] for k in ['bytes','sha256']};paths.append(p)
 inputs={str(p.relative_to(R)):pin(p) for p in paths}
 write(W/'plan.json',dict(commands=commands,outputs=outputs,inputs=inputs,cleanup_pipeline=pipeline,
  proof_scope='Only nine core ordinary OtherRead is removed. Strict/inaccessibleRW remains. First math+DP only; source/no NaN assumptions and no changed raw entry effects.'))
 actions={}
 predecessors={'plugin':None,'math-lower':'plugin','dp-lower':'plugin','dp-cleanup':'dp-lower','math-native':'math-lower','dp-native':'dp-cleanup'}
 for name in commands:
  prior=predecessors[name];extra=[] if prior is None else [str((O/outputs[prior]).relative_to(R))]
  if name.endswith('-lower'):extra.append(str((I/'math.bc').relative_to(R)))
  actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
   limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=[sys.executable,'-B',str(Path(__file__).resolve()),name],
   frozen_inputs=sorted(inputs)+[str((W/'plan.json').relative_to(R))]+extra,
   scope='Bounded local-memory-effect refinement and first DP compiler cleanup, preserving all strict FP/status effects, arithmetic, original raw guards/layout/inputs and work limits. No workload, numerical, speed, provider, frame, scientific or release acceptance.')
 write(W/'actions.json',actions);print(json.dumps(dict(prepared=list(commands),runtime_run=False)))
else:
 plan=doc(W/'plan.json');assert phase in plan['commands']
 for name,expected in plan['inputs'].items():assert pin(R/name)==expected,name
 output=O/plan['outputs'][phase];assert not output.exists()
 predecessors={'plugin':None,'math-lower':'plugin','dp-lower':'plugin','dp-cleanup':'dp-lower','math-native':'math-lower','dp-native':'dp-cleanup'}
 prior=predecessors[phase]
 if prior:
  owner=doc(W/prior/'owner.json');assert owner['passed'] and owner['status']=='terminal'
  product=doc(W/(prior+'-result.json'))['output'];assert pin(R/product['path'])=={k:product[k] for k in ['bytes','sha256']}
 subprocess.run(plan['commands'][phase],check=True)
 result=dict(phase=phase,output=dict(path=str(output.relative_to(R)),**pin(output)),argv=plan['commands'][phase])
 if output.suffix=='.bc':
  ll=output.with_suffix('.ll');subprocess.run([str(D/'llvm-dis'),str(output),'-o',str(ll)],check=True)
  text=ll.read_text();definitions=re.findall(r'^define .*$',text,re.M);attrs=dict(re.findall(r'^attributes #(\d+) = (.*)$',text,re.M))
  assert all('strictfp' in attrs[re.search(r'#(\d+)',line)[1]] for line in definitions)
  assert not re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',text,re.M)
  ordinary=Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b',text,re.M))
  assert not ordinary,ordinary
  core=Counter(re.findall(r'\bcall\s+.*?@('+PREFIX+r'[^ (]+)\(',text));assert sum(core.values())==(0 if phase=='math-lower' else 1581)
  result.update(ir=dict(path=str(ll.relative_to(R)),**pin(ll)),definitions=len(definitions),core_calls=dict(core),
   constrained_calls=dict(Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(',text))),
   memory_instruction_counts=dict(Counter(re.findall(r'^\s+(?:%[^=]+ = )?(load|store|alloca)\b',text,re.M))),
   masked_store_calls=text.count('call void @llvm.masked.store'),ordinary_FP_forbidden=True)
 elif phase.endswith('-native'):
  asm=subprocess.check_output(['objdump','-drw',str(output)],text=True)
  instructions=re.findall(r'^\s*[0-9a-f]+:\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9.]*)',asm,re.M)
  assert instructions and not [x for x in instructions if x.startswith(('vfm','vfnm'))]
  refs=Counter(re.findall(r'R_X86_64_PLT32\s+('+PREFIX+r'\w+)-0x4',asm));assert sum(refs.values())==(0 if phase=='math-native' else 1581)
  result.update(instructions=len(instructions),core_relocations=dict(refs),emitted_fused_instructions=[])
 write(W/(phase+'-result.json'),result);print(json.dumps({k:result[k] for k in ['phase','output','instructions','memory_instruction_counts','masked_store_calls'] if k in result}))
