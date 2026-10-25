"""Replace only the three corrected raw stages and private provider identities."""
from pathlib import Path
import difflib,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider';S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8';V=R/'out/native-cpu-predicated-refinement/provider'
J=lambda p:json.loads(p.read_text());pin=lambda b:dict(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
canonical=lambda v:pin(json.dumps(v,sort_keys=True,separators=(',',':')).encode())
assert not (P/'provider.patch').exists();old=J(V/'source-contract.json')
before={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()};assert canonical(before)==old['full_candidate_source_identity']
compiled=J(W/'ispc-dense-compile-result.json');assert compiled['all_objects_compiled'] and compiled['new_raw_compilations']==3 and compiled['new_CPP_compilations']==0
assert J(W/'dense-code-inspection.json')['passed'] and J(W/'object-code-inspection-dense.json')['objects']
regression=J(W/'fixture/result.json');assert regression['passed'] and regression['old_stale_tag_arithmetic_counterexamples']==12
assert J(W/'fixture/run/owner.json')['passed'] and J(W/'fixture/build-admission.json')['passed']
objects=old['native_objects'][:4]+[x['object'] for x in compiled['objects']]+old['native_objects'][7:];assert len(objects)==8
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
texts={};name='src/gallium/targets/lavapipe/meson.build';text=(S/name).read_text()
for a,b in zip(old['native_objects'][4:7],objects[4:7]):
 assert text.count(str(R/a['path']))==1;text=text.replace(str(R/a['path']),str(R/b['path']),1)
texts[name]=text
for name,a,b in [('src/gallium/drivers/llvmpipe/lp_screen.c','Sirius-native-cpu-predicated-refinement-v35','Sirius-native-cpu-reachable-tag-writers-v36'),('src/gallium/frontends/lavapipe/lvp_device.c','SIR-CPU-25.2.8AK','SIR-CPU-25.2.8AL')]:
 text=(S/name).read_text();assert text.count(a)==1;texts[name]=text.replace(a,b,1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
 original=(S/name).read_text();patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,fromfile='a/'+name,tofile='b/'+name));after[name]=pin(text.encode());modified[name]=dict(basis=before[name],candidate=after[name])
inputs={**old['source_inputs'],**compiled['source_inputs']}
for row in objects:inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
for p in [V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',W/'ispc-derivation.json',W/'ispc-dense-compile-result.json',W/'dense-code-inspection.json',W/'object-code-inspection-dense.json',W/'compile/ispc-dense/owner.json',W/'writer-fixture-source.json',W/'fixture/result.json',W/'fixture/run/owner.json',W/'fixture/build/owner.json',W/'fixture/build-admission.json',W/'fixture/linked-code-admission.json',Path(__file__).resolve()]:inputs[str(p.relative_to(R))]=pin(p.read_bytes())
for name,v in inputs.items():assert pin((R/name).read_bytes())==v,name
contract=dict(basis_source_identity=old['full_candidate_source_identity'],full_candidate_source_identity=canonical(after),modified_files=modified,patch=pin(patch.encode()),
 native_objects=objects,source_inputs=inputs,new_native_compilations=3,new_ISPC_compilations=3,new_CPP_compilations=0,reused_raw_objects=0,
 guard_shared_byte_authority='out/native-cpu-active-prefix/native_packet.hpp',
 superseded_claim='AJ/AK raw writer tag updates were after unconditional returns. Prior coherence claims are refuted by12 original-body stale-tag arithmetic counterexamples. Historical finite outcomes stay scoped observations; do not transfer coherent-tag admission.',
 scope='AL moves only twelve metadata updates before original writer returns. Complete payload stores precede normalization-tag refresh; every partial write invalidates under its original mask/index. Exact source inverse restores AJ; all9380ordered tuples/math/public layouts/eligibility/FP/ISA/requests/fields/RHS/trials/quartics/count/cache/governors/order/rollback remain. Exact AK arithmetic/adapter and scalar3/pure admission reused. Actual-body regression passes190512writer calls/22861440bank/4572288tag words and12counterexamples against original scalar. Original build-reader failure remains explicitly admitted without compiler/link retry. Current provider/code/full38/C8/affected controls still required; no speed/science/release acceptance.')
P.mkdir(exist_ok=True);(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps(dict(modified_files=3,new_raw_objects=3,CPP_math_recompiled=False,shared_Mesa_modified=False)))
