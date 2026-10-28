"""Read retained successful lowering after a pre-inlining count receipt defect."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,subprocess,sys

R=Path.cwd();W=Path(__file__).resolve().parent;V=R/'out/native-cpu-strict-ir-linking';B=R/'bin/linux-gcc/native-cpu-canonical-closed-dp'
P=B/'dp.value.bc';EXPECTED=dict(bytes=10646600,sha256='3f69125b95e0a6282a3e7dfa2d7617a0ea933af22dd2aa411314a8b1c76a2cf0')


def pin(p):
 with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())


def doc(p):return json.loads(p.read_text())


def write(p,v):
 with p.open('x') as s:json.dump(v,s,indent=2);s.write('\n')


phase,=sys.argv[1:]
if phase=='prepare':
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 owner=doc(W/'value/owner.json')
 assert owner['status']=='terminal' and not owner['passed'] and owner['child_exit']==1
 assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
 assert pin(P)==EXPECTED
 log=(W/'value/stderr.log').read_text()
 assert 'shared_core_declarations=9 adapters=5 raw_calls=1576 tags=1575' in log
 assert 'FP_inventory_unchanged=1' in log and "sum(result['value_calls'].values())" in log
 paths=[Path(__file__).resolve(),W/'pipeline.py',W/'value/owner.json',W/'value/stderr.log',W/'scope-result.json',V/'closed_dp_value_contract.json',V/'closed_dp_scope_contract.json',P,B/'dp.value.ll']
 inputs={str(p.relative_to(R)):pin(p) for p in paths}
 write(W/'value-audit-plan.json',dict(inputs=inputs,original_failed_owner_preserved=True,lowering_repeated=False))
 actions=doc(W/'actions.json');assert 'value-audit' not in actions
 (W/'value-audit-prior-actions.json').write_bytes((W/'actions.json').read_bytes())
 actions['value-audit']=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
  limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=[sys.executable,'-B',str(Path(__file__).resolve()),'audit'],
  frozen_inputs=sorted(inputs)+[str((W/'value-audit-plan.json').relative_to(R))],
  scope='Read-only validation of retained verified lowering after Python counted typed calls before adapter inlining. Original failure intact; no conversion, optimization, native/workload or qualification rerun.')
 (W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
 print(json.dumps(dict(prepared='value-audit',lowering_repeated=False)))
else:
 assert phase=='audit'
 for name,expected in doc(W/'value-audit-plan.json')['inputs'].items():assert pin(R/name)==expected,name
 assert pin(P)==EXPECTED
 text=(B/'dp.value.ll').read_text();contract=doc(V/'closed_dp_value_contract.json')
 expected=contract['exact_calls_by_recipe']
 pointer=Counter(re.findall(r'\bcall\s+.*?@(sirius_closed_scope_v4_[^ (]+)\(',text))
 typed=Counter(re.findall(r'\bcall\s+.*?@(sirius_local_value_[^ (]+)\(',text))
 assert pointer==Counter({'sirius_closed_scope_v4_'+k:v for k,v in expected.items()})
 assert typed==Counter({'sirius_local_value_'+k:1 for k in expected})
 attrs=dict(re.findall(r'^attributes #(\d+) = (.*)$',text,re.M))
 adapters=re.findall(r'^define internal .*@sirius_closed_scope_v4_.*$',text,re.M)
 assert len(adapters)==5 and all('alwaysinline' in attrs[re.search(r'#(\d+)',h)[1]] and 'strictfp' in attrs[re.search(r'#(\d+)',h)[1]] for h in adapters)
 constrained=Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(',text))
 assert constrained==Counter(doc(W/'scope-result.json')['constrained_calls'])
 assert not re.findall(r'^\s+%[^=]+ = (?:fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b',text,re.M)
 result=dict(phase='value',output=dict(path=str(P.relative_to(R)),**pin(P)),input=contract['input'],scope_contract=pin(V/'closed_dp_scope_contract.json'),
  ir=dict(path=str((B/'dp.value.ll').relative_to(R)),**pin(B/'dp.value.ll')),pointer_calls=dict(pointer),value_calls=dict(typed),constrained_calls=dict(constrained),
  source_value_contract=pin(V/'closed_dp_value_contract.json'),original_failed_owner=pin(W/'value/owner.json'),
  accepted_artifact_owner=str((W/'value-audit/owner.json').relative_to(R)),lowering_repeated=False,
  scope='Typed5 wrapper calls and1576 original pointer calls coexist before always-inlining. Existing opt already verified the full newly pinned transformation and FP inventory. This reader preserves the failed receipt and verifies the artifact; no runtime/numerical/speed qualification.')
 write(W/'value-result.json',result);print(json.dumps(dict(existing_lowering_verified=True,original_pointer_calls=sum(pointer.values()),preinline_typed_adapter_calls=sum(typed.values()),lowering_repeated=False)))
