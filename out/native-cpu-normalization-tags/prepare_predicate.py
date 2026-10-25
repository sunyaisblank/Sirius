"""Reuse the compiled driver/math; link and run only the shared predicate check."""
from pathlib import Path
import hashlib,json,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'predicate';B=R/'bin/linux-gcc/native-cpu-normalization-tags/predicate'
J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
mode,=sys.argv[1:];assert mode in ('link','run') and not subprocess.check_output(['git','status','--porcelain']).strip()
c=J(W/'cpp-build-result.json');i=J(W/'ispc-dense-compile-result.json');assert c['passed'] and i['all_objects_compiled'] and i['new_compilations']==4
for p in (W/'cpp/build/owner.json',W/'cpp/run/owner.json',W/'compile/ispc-dense/owner.json'):
 o=J(p);assert o['passed'] and o['status']=='terminal' and o['source_unchanged'] and o['whole_inputs_unchanged'] and not o['remaining_owned_processes'] and not o['cleanup_errors']
assert J(W/'cpp/run/stdout.log')['compared_words']==44452800
objects=[x['object'] for x in c['compilations']]+[i['predicate_compilation']['object']]
assert len(objects)==3 and all(pin(R/x['path'])=={k:x[k] for k in ('bytes','sha256')} for x in objects)
inputs={**J(W/'prepared-cpp.json')['inputs'],**i['source_inputs']}
for x in objects:inputs[x['path']]={k:x[k] for k in ('bytes','sha256')}
for p in (W/'cpp-build-result.json',W/'ispc-dense-compile-result.json',W/'cpp/run/owner.json',W/'compile/ispc-dense/owner.json',W/'prepare_predicate.py',W/'predicate_task.py',P/'run_owned.py',P/'owned_guard.py',R/'out/native-cpu-field-compatibility/elf_proof.py'):
 inputs[str(p.relative_to(R))]=pin(p)
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
if mode=='link':
 assert not B.exists() and not (P/'prepared.json').exists() and not (P/'actions.json').exists()
 argv=[J(W/'prepared-cpp.json')['link'][0],*[str(R/x['path']) for x in objects],'-o',str(B/'abi-check')]
 prepared={'revision':revision,'inputs':inputs,'objects':objects,'argv':argv,'executable':str((B/'abi-check').relative_to(R)),'scope':'One fixture-only link reuses exact already-compiled CPP driver/math and shared-predicate ISPC object. No new C++/ISPC/math compilation. Predicate-only mode exercises the exact shared ISPC integer predicate against the independent std::nextafter scalar predicate; no arithmetic workload repeats or performance/qualification claim.'}
 (P/'prepared.json').write_text(json.dumps(prepared,indent=2)+'\n');actions={};kind='build_only'
else:
 owner=J(P/'link/owner.json');assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
 built=J(P/'link-result.json');assert built['passed'] and built['new_links']==1 and built['new_compilations']==0
 inputs.update(built['runtime_libraries']);inputs[built['executable']['path']]={k:built['executable'][k] for k in ('bytes','sha256')}
 for p in (P/'prepared.json',P/'link-result.json',P/'link/owner.json'):inputs[str(p.relative_to(R))]=pin(p)
 actions=J(P/'actions.json');assert 'run' not in actions and not (P/'run').exists();kind='diagnostic'
for n,v in inputs.items():assert pin(R/n)==v,n
inputs[str((P/'prepared.json').relative_to(R))]=pin(P/'prepared.json')
a={'revision':revision,'kind':kind,'expected_exit':0,'software_provider':False,'limits':{'seconds':300,'rss_kib':4194304,'reserve_kib':2097152},'frozen_inputs':sorted(inputs),'argv':['/usr/bin/python3','-B',str(W/'predicate_task.py'),mode], 'scope':'Shared normalization-predicate fixture only. Exact original49tuples, broadcast/mixed packets, masks, destination columns, inactive tags and unchanged input words, independent std::nextafter reference. No arithmetic test retry or performance/science/release qualification.'}
actions[mode]=a;(P/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps({'prepared':mode,'new_compilations':0,'executed':False}))
