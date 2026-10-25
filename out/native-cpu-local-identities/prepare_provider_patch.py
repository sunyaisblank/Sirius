"""Link only three new raw objects, retaining every AM C++ math byte."""
from pathlib import Path
import difflib,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider';S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8';V=R/'out/native-cpu-value-fact-stores/provider'
J=lambda p:json.loads(p.read_text());pin=lambda b:dict(bytes=len(b),sha256=hashlib.sha256(b).hexdigest());canonical=lambda v:pin(json.dumps(v,sort_keys=True,separators=(',',':')).encode())
assert not (P/'provider.patch').exists();old=J(V/'source-contract.json')
before={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()};assert canonical(before)==old['full_candidate_source_identity']
compiled=J(W/'ispc-dense-compile-result.json');assert compiled['all_objects_compiled'] and compiled['new_raw_compilations']==3 and compiled['new_CPP_compilations']==0
assert J(W/'dense-code-inspection.json')['passed'] and J(V.parent/'cpp-code-inspection.json')['passed']
reg=J(W/'fixture/result.json');assert reg['passed'] and reg['calls']==555660 and reg['changed_operator_calls']==333396
assert J(W/'fixture/run/owner.json')['passed'] and J(W/'fixture/build/owner.json')['passed']
objects=old['native_objects'][:4]+[x['object'] for x in compiled['objects']]+[compiled['merged_adapter']];assert len(objects)==8 and objects[-1]==old['native_objects'][-1]
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
texts={};name='src/gallium/targets/lavapipe/meson.build';text=(S/name).read_text()
for a,b in zip(old['native_objects'][4:7],objects[4:7]):
 assert text.count(str(R/a['path']))==1;text=text.replace(str(R/a['path']),str(R/b['path']),1)
texts[name]=text
for name,a,b in [('src/gallium/drivers/llvmpipe/lp_screen.c','Sirius-native-cpu-value-fact-stores-v37','Sirius-native-cpu-local-identities-v38'),('src/gallium/frontends/lavapipe/lvp_device.c','SIR-CPU-25.2.8AM','SIR-CPU-25.2.8AN')]:
 text=(S/name).read_text();assert text.count(a)==1;texts[name]=text.replace(a,b,1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
 original=(S/name).read_text();patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,fromfile='a/'+name,tofile='b/'+name));after[name]=pin(text.encode());modified[name]=dict(basis=before[name],candidate=after[name])
inputs={**old['source_inputs'],**compiled['source_inputs']}
for row in objects:inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
for p in [V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',W/'ispc-derivation.json',W/'ispc-dense-compile-result.json',W/'dense-code-inspection.json',W/'object-code-inspection-dense.json',W/'compile/ispc-dense/owner.json',W/'fixture/result.json',W/'fixture/run/owner.json',W/'fixture/build/owner.json',W/'fixture/build-result.json',W/'fixture/linked-code.json',W/'identity_store.isph',Path(__file__).resolve()]:inputs[str(p.relative_to(R))]=pin(p.read_bytes())
for name,v in inputs.items():assert pin((R/name).read_bytes())==v,name
contract=dict(basis_source_identity=old['full_candidate_source_identity'],full_candidate_source_identity=canonical(after),modified_files=modified,patch=pin(patch.encode()),native_objects=objects,source_inputs=inputs,new_native_compilations=3,new_ISPC_compilations=3,new_CPP_compilations=0,reused_raw_objects=0,guard_shared_byte_authority='out/native-cpu-active-prefix/native_packet.hpp',scope='AN moves only original AM complete-fact identity/refusal stores before the same C++ ABI calls in three raw stages. Integer payload copies/sign bits, original precedence/valid/signed-zero rules and alias snapshots remain; every unknown/legacy/general active packet takes unchanged AM helper. Same metadata arrays/writers, exact9380tuples and all general expressions/refinements/guards/FP/ISA/count/cache/governors/fields/RHS/trials/quartics/scheduler/rollback. All five other native objects and host remain exact. One actual common-include fixture passes555660scalar-reference calls/333396changed-operator calls and fullbank/metadata comparisons. Actual provider/full38/C8 and any adoption benefit/control/frame gates remain separate; no current performance/science/release claim.')
(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()};print(json.dumps(dict(modified_files=3,new_runtime_raw_objects=3,new_runtime_CPP_helper=0,shared_Mesa_modified=False)))
