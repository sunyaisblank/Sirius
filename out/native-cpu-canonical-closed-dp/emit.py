"""Derive and emit the first closed private DP source under original limits."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,subprocess,sys

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking'
B=R/'bin/linux-gcc/native-cpu-canonical-closed-dp'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def write(p,v):
 with p.open('x') as s:json.dump(v,s,indent=2);s.write('\n')


phase,=sys.argv[1:];contract=json.loads((V/'canonical_private_dp_contract.json').read_text())
if phase=='prepare':
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 assert not B.exists()
 paths=[*W.glob('*.py'),V/'canonical_private_dp.py',V/'canonical_private_dp_contract.json',
  R/'out/native-cpu-row-ownership/cold-frame/run_owned.py',R/'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
 inputs={str(p.relative_to(R)):pin(p) for p in paths}|contract['source_inputs']
 for name,expected in inputs.items():assert pin(R/name)==expected,name
 commands=dict(derive=[sys.executable,'-B',str(V/'canonical_private_dp.py'),'derive'],emit=contract['prospective_emit_command'])
 write(W/'plan.json',dict(commands=commands,inputs=inputs))
 actions={}
 for name in commands:
  extra=[] if name=='derive' else [str((B/'ispc-source/dopri_phase.ispc').relative_to(R))]
  actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
   limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=[sys.executable,'-B',str(Path(__file__).resolve()),name],
   frozen_inputs=sorted(inputs)+[str((W/'plan.json').relative_to(R))]+extra,
   scope='First closed private canonical DP source/compiler emission; positive immutable-program job/full FP scope required. Same original tuples/layout/allocation/guards/refusals/arithmetic and work limits. No strict conversion, executable, numerical, speed, provider/frame/scientific/release qualification.')
 write(W/'actions.json',actions);print(json.dumps(dict(prepared=list(commands),runtime_run=False)))
else:
 plan=json.loads((W/'plan.json').read_text());assert phase in plan['commands']
 for name,expected in plan['inputs'].items():assert pin(R/name)==expected,name
 if phase=='emit':
  owner=json.loads((W/'derive/owner.json').read_text());assert owner['passed'] and owner['status']=='terminal'
  assert pin(B/'ispc-source/dopri_phase.ispc')==contract['derived_source']
  assert not (B/'dopri_phase.ll').exists()
 subprocess.run(plan['commands'][phase],check=True)
 if phase=='derive':
  write(W/'derive-result.json',dict(output=dict(path=str((B/'ispc-source/dopri_phase.ispc').relative_to(R)),**pin(B/'ispc-source/dopri_phase.ispc')),source_proof=contract['source_proof']))
 else:
  p=B/'dopri_phase.ll';text=p.read_text();definitions=re.findall(r'^define .*$',text,re.M)
  ordinary=Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b',text,re.M))
  calls=Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(',text));predicates=Counter(re.findall(r'= fcmp (\w+) ',text))
  fast=re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',text,re.M);assert not fast
  result=dict(output=dict(path=str(p.relative_to(R)),**pin(p)),definitions=len(definitions),ordinary_fp_operations=dict(ordinary),fp_instruction_flags={},
   comparison_predicates=dict(predicates),call_inventory=dict(calls),math_calls={k:v for k,v in calls.items() if k.startswith('sirius_native_v4_')},
   definition_headers=definitions,scope='First ordinary ISPC IR only. Source canonical guard/order proof retained; no strict/body/native/adoption qualification.')
  write(W/'inventory.json',result);print(json.dumps({k:result[k] for k in ['output','definitions','ordinary_fp_operations','comparison_predicates','math_calls']}))
