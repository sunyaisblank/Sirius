"""Compile exact derived strict/value passes; never relax immutable input pins."""
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


phase,=sys.argv[1:]
mode=phase.removeprefix('prepare-').removesuffix('-plugin')
assert mode in ['scope','value']
contract_path=V/('closed_dp_scope_contract.json' if mode=='scope' else 'closed_dp_value_contract.json')
contract=doc(contract_path)
if phase.startswith('prepare-'):
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 source=R/contract['output']['path'];assert pin(source)=={k:contract['output'][k] for k in ['bytes','sha256']}
 config=shlex.split(subprocess.check_output([str(D/'llvm-config'),'--cxxflags'],text=True))
 outputs={mode+'-plugin':mode+'.so',mode:'dp.'+('strict' if mode=='scope' else 'value')+'.bc'}
 commands={mode+'-plugin':[str(D/'clang++'),*config,'-O2','-fPIC','-shared',str(source),str(SDK/'sha256-only.a'),'-Wl,--exclude-libs,ALL','-o',str(B/outputs[mode+'-plugin'])],
  mode:[str(D/'opt'),'-load-pass-plugin='+str(B/outputs[mode+'-plugin']),'-passes='+contract['pass_name']+',verify','-verify-each',str(R/contract['input']['path']),'-o',str(B/outputs[mode])]}
 paths=[*W.glob('*.py'),V/'derive_closed_dp_passes.py',source,contract_path,SDK/'sha256-only.a',V/'toolchain.json']
 inputs={str(p.relative_to(R)):pin(p) for p in paths}|contract['source_inputs']
 for n,p in inputs.items():assert pin(R/n)==p,n
 for name in ['clang++','opt','llvm-config','llvm-dis']:
  tool=doc(V/'toolchain.json')['tools'][name];p=R/tool['resolved'];assert pin(p)=={k:tool[k] for k in ['bytes','sha256']};inputs[str(p.relative_to(R))]=pin(p)
 plan_path=W/(mode+'-plan.json');write(plan_path,dict(commands=commands,outputs=outputs,inputs=inputs,contract=pin(contract_path)))
 actions=doc(W/'actions.json')
 # Retain exact historical action bytes before appending new compiler actions.
 previous=W/(mode+'-prior-actions.json');previous.write_bytes((W/'actions.json').read_bytes())
 for name in commands:
  assert name not in actions
  extra=[] if name.endswith('-plugin') else [str((B/outputs[mode+'-plugin']).relative_to(R))]
  actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
   limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=[sys.executable,'-B',str(Path(__file__).resolve()),name],
   frozen_inputs=sorted(inputs)+[str(plan_path.relative_to(R))]+extra,
   scope='Exact separately named closed-DP strict/value implementation. Original pinned conversion/bank proof reused only for newly bound actual source with same initialization/admitted job/FP scope. All types/rounding/strict effects/guards preserved; no force/trust path, new math body, runtime/speed/provider/frame/scientific/release claim.')
 (W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps(dict(prepared=list(commands),runtime_run=False)))
else:
 plan=doc(W/(mode+'-plan.json'));assert phase in plan['commands']
 for n,p in plan['inputs'].items():assert pin(R/n)==p,n
 if phase==mode:
  o=doc(W/(mode+'-plugin')/'owner.json');assert o['passed'] and o['status']=='terminal'
  product=doc(W/(mode+'-plugin-result.json'))['output'];assert pin(R/product['path'])=={k:product[k] for k in ['bytes','sha256']}
 output=B/plan['outputs'][phase];assert not output.exists()
 subprocess.run(plan['commands'][phase],check=True)
 result=dict(phase=phase,output=dict(path=str(output.relative_to(R)),**pin(output)),input=contract['input'],
  scope_contract=pin(V/'closed_dp_scope_contract.json'),argv=plan['commands'][phase])
 if phase==mode:
  ll=output.with_suffix('.ll');subprocess.run([str(D/'llvm-dis'),str(output),'-o',str(ll)],check=True)
  text=ll.read_text();headers=re.findall(r'^define .*$',text,re.M);attrs=dict(re.findall(r'^attributes #(\d+) = (.*)$',text,re.M))
  assert headers and all('strictfp' in attrs[re.search(r'#(\d+)',h)[1]] for h in headers)
  ordinary=Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b',text,re.M));assert not ordinary
  assert not re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',text,re.M)
  result.update(ir=dict(path=str(ll.relative_to(R)),**pin(ll)),definitions=len(headers),
   constrained_calls=dict(Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(',text))),
   value_calls=dict(Counter(re.findall(r'\bcall\s+.*?@(sirius_local_value_[^ (]+)\(',text))))
  if mode=='value':assert sum(result['value_calls'].values())==sum(contract['exact_calls_by_recipe'].values())
 write(W/(phase+'-result.json'),result);print(json.dumps({k:result[k] for k in ['phase','output','definitions','value_calls'] if k in result}))
