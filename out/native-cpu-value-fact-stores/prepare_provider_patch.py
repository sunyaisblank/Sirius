"""Replace runtime raw3/helper with complete-fact identity stores only."""
from pathlib import Path
import difflib,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider';S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8';V=R/'out/native-cpu-reachable-tag-writers/provider'
J=lambda p:json.loads(p.read_text());pin=lambda b:dict(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
canonical=lambda v:pin(json.dumps(v,sort_keys=True,separators=(',',':')).encode())
assert not (P/'provider.patch').exists();old=J(V/'source-contract.json')
before={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()};assert canonical(before)==old['full_candidate_source_identity']
compiled=J(W/'ispc-dense-compile-result.json');cpp=J(W/'cpp-build-result.json')
assert compiled['all_objects_compiled'] and compiled['new_raw_compilations']==3 and cpp['passed'] and cpp['new_relocatable_links']==1
assert J(W/'dense-code-inspection.json')['passed'] and J(W/'cpp-code-inspection.json')['passed']
math=J(W/'cpp/run/stdout.log');assert math['passed'] and math['calls']==555660 and math['compared_words']==66679200 and math['tag_words']==13335840
regression=J(W/'fixture/result.json');assert regression['passed'] and regression['fixed_writer_calls']==254016 and regression['old_stale_tag_arithmetic_counterexamples']==24
assert J(W/'fixture/run/owner.json')['passed'] and J(W/'fixture/build/owner.json')['passed']
objects=old['native_objects'][:4]+[x['object'] for x in compiled['objects']]+[cpp['merged_adapter']];assert len(objects)==8
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
texts={};name='src/gallium/targets/lavapipe/meson.build';text=(S/name).read_text()
for a,b in zip(old['native_objects'][4:],objects[4:]):
 assert text.count(str(R/a['path']))==1;text=text.replace(str(R/a['path']),str(R/b['path']),1)
texts[name]=text
for name,a,b in [('src/gallium/drivers/llvmpipe/lp_screen.c','Sirius-native-cpu-reachable-tag-writers-v36','Sirius-native-cpu-value-fact-stores-v37'),('src/gallium/frontends/lavapipe/lvp_device.c','SIR-CPU-25.2.8AL','SIR-CPU-25.2.8AM')]:
 text=(S/name).read_text();assert text.count(a)==1;texts[name]=text.replace(a,b,1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
 original=(S/name).read_text();patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,fromfile='a/'+name,tofile='b/'+name));after[name]=pin(text.encode());modified[name]=dict(basis=before[name],candidate=after[name])
inputs={**old['source_inputs'],**compiled['source_inputs'],**J(W/'prepared-cpp.json')['inputs']}
for row in objects:inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
for p in [V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',W/'ispc-derivation.json',W/'ispc-dense-compile-result.json',W/'dense-code-inspection.json',W/'object-code-inspection-dense.json',W/'compile/ispc-dense/owner.json',W/'writer-fixture-source.json',W/'fixture/result.json',W/'fixture/run/owner.json',W/'fixture/build/owner.json',W/'fixture/build-result.json',W/'fixture/linked-code.json',W/'cpp-source-proof.json',W/'cpp-build-result.json',W/'cpp-code-inspection.json',W/'cpp/run/owner.json',W/'cpp/run/stdout.log',W/'value_facts.isph',Path(__file__).resolve()]:inputs[str(p.relative_to(R))]=pin(p.read_bytes())
for name,v in inputs.items():assert pin((R/name).read_bytes())==v,name
contract=dict(basis_source_identity=old['full_candidate_source_identity'],full_candidate_source_identity=canonical(after),modified_files=modified,patch=pin(patch.encode()),native_objects=objects,source_inputs=inputs,
 new_native_compilations=4,new_ISPC_compilations=3,new_CPP_compilations=1,reused_raw_objects=0,guard_shared_byte_authority='out/native-cpu-active-prefix/native_packet.hpp',
 scope='AM reuses exact zero/one/minus-one/nonzero-valid facts in the existing private metadata word. Complete words alone authorize original identity/refusal bank stores, including zero precedence/copy signs/valid2; selected operands are snapshotted before aliasing writes, zero/refusal packets avoid payload loads. Unknown/legacy metadata falls back conservatively. General expressions/order are extracted exactly from AK; original untagged/Divide/Sqrt and full fallback remain. All partial writes invalidate before returns, complete writes refresh all exact facts. Same public5words/context arrays/canonical9380tuples/guards/precision/FP/ISA/count/cache/governors/fields/RHS/trials/quartics/order/rollback/scheduler. Fresh555660call/66679200bank/13335840metadata scalar comparison and254016actual-writer calls/24legacy-rich stale-tag witnesses pass. Inherited writer fixture scope text says AK math, but its sealed actual helper identity and full code join are the new AM object; old records are not rewritten. Actual provider/full38/C8/original performance/control/frame gates remain required; no speed/science/release claim.')
P.mkdir(exist_ok=True);(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps(dict(modified_files=3,new_runtime_raw_objects=3,new_runtime_CPP_helper=1,shared_Mesa_modified=False)))
