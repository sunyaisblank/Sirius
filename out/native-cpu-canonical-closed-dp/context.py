"""One whole-owned-context normalization, cleanup and native compiler gate."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,shlex,subprocess,sys

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking'
B=R/'bin/linux-gcc/native-cpu-canonical-closed-dp';SDK=R/'bin/linux-gcc/native-cpu-strict-ir-linking';D=SDK/'toolchain/bin'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def doc(p):return json.loads(p.read_text())


def write(p,v):
 with p.open('x') as s:json.dump(v,s,indent=2);s.write('\n')


phase,=sys.argv[1:];contract=doc(V/'closed_dp_context_contract.json')
if phase=='prepare':
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 source=R/contract['output']['path'];assert pin(source)=={k:contract['output'][k] for k in ['bytes','sha256']}
 config=shlex.split(subprocess.check_output([str(D/'llvm-config'),'--cxxflags'],text=True))
 pipeline='function(sroa,mem2reg,instcombine,simplifycfg,early-cse<memssa>,dse),globaldce,verify'
 commands={'context-plugin':[str(D/'clang++'),*config,'-O2','-fPIC','-shared',str(source),str(SDK/'sha256-only.a'),'-Wl,--exclude-libs,ALL','-o',str(B/'context.so')],
  'context':[str(D/'opt'),'-load-pass-plugin='+str(B/'context.so'),'-passes='+contract['pass_name']+',verify','-verify-each',str(R/contract['input']['path']),'-o',str(B/'dp.context.bc')],
  'cleanup':[str(D/'opt'),'-passes='+pipeline,'-verify-each','-mcpu=skylake-avx512',str(B/'dp.context.bc'),'-o',str(B/'dp.closed.bc')],
  'native':[str(D/'llc'),'-O2','-filetype=obj','-relocation-model=pic','-mcpu=skylake-avx512',str(B/'dp.closed.bc'),'-o',str(B/'dp.closed.o')]}
 outputs={'context-plugin':'context.so','context':'dp.context.bc','cleanup':'dp.closed.bc','native':'dp.closed.o'}
 paths=[*W.glob('*.py'),source,V/'closed_dp_context_contract.json',V/'derive_closed_dp_context.py',V/'toolchain.json',SDK/'sha256-only.a']
 inputs={str(p.relative_to(R)):pin(p) for p in paths}|contract['source_inputs']
 for n,p in inputs.items():assert pin(R/n)==p,n
 for name in ['clang++','opt','llvm-config','llvm-dis','llc']:
  tool=doc(V/'toolchain.json')['tools'][name];p=R/tool['resolved'];assert pin(p)=={k:tool[k] for k in ['bytes','sha256']};inputs[str(p.relative_to(R))]=pin(p)
 write(W/'context-plan.json',dict(commands=commands,outputs=outputs,inputs=inputs,cleanup_pipeline=pipeline))
 actions=doc(W/'actions.json');(W/'context-prior-actions.json').write_bytes((W/'actions.json').read_bytes())
 prior={'context-plugin':None,'context':'context-plugin','cleanup':'context','native':'cleanup'}
 for name in commands:
  assert name not in actions
  extra=[] if not prior[name] else [str((B/outputs[prior[name]]).relative_to(R))]
  actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
   limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=[sys.executable,'-B',str(Path(__file__).resolve()),name],
   frozen_inputs=sorted(inputs)+[str((W/'context-plan.json').relative_to(R))]+extra,
   scope='Single exact private DP whole-context compiler gate. Inline raw helpers only, keep existing external strict math, normalize only dominated fully initialized120832B local context; no global/other-bank rewrites or scene/arithmetic/precision/work-limit changes. No runtime/provider/benefit/frame/scientific/release acceptance.')
 (W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps(dict(prepared=list(commands),runtime_run=False)))
else:
 plan=doc(W/'context-plan.json');assert phase in plan['commands']
 for n,p in plan['inputs'].items():assert pin(R/n)==p,n
 prior={'context-plugin':None,'context':'context-plugin','cleanup':'context','native':'cleanup'}[phase]
 if prior:
  owner=doc(W/prior/'owner.json');assert owner['passed'] and owner['status']=='terminal'
  product=doc(W/(prior+'-result.json'))['output'];assert pin(R/product['path'])=={k:product[k] for k in ['bytes','sha256']}
 output=B/plan['outputs'][phase];assert not output.exists()
 subprocess.run(plan['commands'][phase],check=True)
 result=dict(phase=phase,output=dict(path=str(output.relative_to(R)),**pin(output)),argv=plan['commands'][phase])
 if output.suffix=='.bc':
  ll=output.with_suffix('.ll');subprocess.run([str(D/'llvm-dis'),str(output),'-o',str(ll)],check=True)
  text=ll.read_text();headers=re.findall(r'^define .*$',text,re.M);attrs=dict(re.findall(r'^attributes #(\d+) = (.*)$',text,re.M))
  assert len(headers)==1 and '@sirius_closed_dp_packet(' in headers[0]
  assert 'strictfp' in attrs[re.search(r'#(\d+)',headers[0])[1]]
  assert not re.findall(r'^\s+%[^=]+ = (?:fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b',text,re.M)
  assert not re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',text,re.M)
  calls=Counter(re.findall(r'\bcall\s+.*?@(sirius_local_value_[^ (]+)\(',text))
  tagged={k:v for k,v in calls.items() if k.startswith('sirius_local_value_tagged_')}
  assert tagged=={'sirius_local_value_tagged_add8':585,'sirius_local_value_tagged_subtract8':284,'sirius_local_value_tagged_multiply8':661,'sirius_local_value_tagged_divide8':45}
  context_allocas=re.findall(r'^\s+%[^=]+ = alloca %v8_varying_KernelContext_0\b',text,re.M)
  result.update(ir=dict(path=str(ll.relative_to(R)),**pin(ll)),definitions=len(headers),value_calls=dict(calls),
   constrained_calls=dict(Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(',text))),owned_context_allocas=len(context_allocas),
   memory_instruction_counts=dict(Counter(re.findall(r'^\s+(?:%[^=]+ = )?(load|store|alloca)\b',text,re.M))),masked_store_calls=text.count('call void @llvm.masked.store'))
  if phase=='cleanup':assert len(context_allocas)==0,'Required context promotion gate failed; stop before native/runtime.'
 elif phase=='native':
  asm=subprocess.check_output(['objdump','-drw',str(output)],text=True)
  instructions=re.findall(r'^\s*[0-9a-f]+:\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9.]*)',asm,re.M)
  assert instructions and not [x for x in instructions if x.startswith(('vfm','vfnm'))]
  result.update(instructions=len(instructions),emitted_fused_instructions=[],core_relocations=dict(Counter(re.findall(r'R_X86_64_PLT32\s+(sirius_local_value_\w+)-0x4',asm))))
 write(W/(phase+'-result.json'),result);print(json.dumps({k:result[k] for k in ['phase','output','owned_context_allocas','instructions','memory_instruction_counts'] if k in result}))
