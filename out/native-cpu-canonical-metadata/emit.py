"""One source derivation and first DP emission, before any coupled qualification."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,subprocess,sys

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking'
B=R/'bin/linux-gcc/native-cpu-canonical-metadata'


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def write(p,value):
 with p.open('x') as s:json.dump(value,s,indent=2);s.write('\n')


phase,=sys.argv[1:]
contract=json.loads((V/'canonical_metadata_contract.json').read_text())
if phase=='prepare':
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 assert not B.exists()
 paths=[*W.glob('*.py'),V/'canonical_metadata.py',V/'canonical_metadata_contract.json',
        R/'out/native-cpu-row-ownership/cold-frame/run_owned.py',R/'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
 inputs={str(p.relative_to(R)):pin(p) for p in paths}|contract['source_inputs']
 for n,p in inputs.items():assert pin(R/n)==p,n
 write(W/'plan.json',dict(inputs=inputs,commands=dict(derive=[sys.executable,'-B',str(V/'canonical_metadata.py'),'derive'],emit=contract['prospective_emit_commands']['dopri_phase']),
  scope='Pinned source insertion/macro change and original ISPC emit flags. First DP structural gate before spending on Transport/Endpoint or any runtime.'))
 actions={}
 for name in ['derive','emit']:
  extra=[] if name=='derive' else [str((B/'ispc-source'/f'{s}.ispc').relative_to(R)) for s in ['dopri_phase','transport','endpoint']]
  actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
   limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=[sys.executable,'-B',str(Path(__file__).resolve()),name],
   frozen_inputs=sorted(inputs)+[str((W/'plan.json').relative_to(R))]+extra,
   scope='First canonical metadata source/DP compiler gate; retain exact original tables/tuples/fallback/arithmetic and original work limits. No numerical, strict-conversion, provider, speed, frame, scientific or release acceptance.')
 write(W/'actions.json',actions)
 print(json.dumps(dict(prepared=['derive','emit'],runtime_run=False)))
else:
 plan=json.loads((W/'plan.json').read_text())
 for n,p in plan['inputs'].items():assert pin(R/n)==p,n
 assert phase in plan['commands']
 if phase=='emit':
  owner=json.loads((W/'derive/owner.json').read_text());assert owner['passed'] and owner['status']=='terminal'
  for stage,expected in contract['derived_sources'].items():assert pin(B/'ispc-source'/(stage+'.ispc'))==expected
  assert not (B/'dopri_phase.ll').exists()
 subprocess.run(plan['commands'][phase],check=True)
 if phase=='derive':
  write(W/'derive-result.json',dict(generated={stage:dict(path=str((B/'ispc-source'/(stage+'.ispc')).relative_to(R)),**pin(B/'ispc-source'/(stage+'.ispc'))) for stage in contract['derived_sources']},source_inverse_exact=True,compiled=False))
 else:
  output=B/'dopri_phase.ll';text=output.read_text();definitions=re.findall(r'^define .*$',text,re.M)
  ordinary=Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b',text,re.M))
  calls=Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(',text))
  fast=re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',text,re.M)
  assert not fast
  write(W/'emit-result.json',dict(output=dict(path=str(output.relative_to(R)),**pin(output)),definitions=len(definitions),ordinary_FP=dict(ordinary),call_inventory=dict(calls),
   scope='Original ISPC ordinary-FP IR only, not yet strict-converted or executable. Actual fast-path metadata/root/context structure still must be checked; no numerical, native or speed acceptance.'))
  print(json.dumps(dict(emitted=str(output.relative_to(R)),definitions=len(definitions),ordinary_FP=dict(ordinary))))
