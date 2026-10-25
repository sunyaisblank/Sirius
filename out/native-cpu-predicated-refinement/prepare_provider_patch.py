"""Stage only the changed helper object and private provider identities."""
from pathlib import Path
import difflib,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider'
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8';V=R/'out/native-cpu-normalization-tags/provider'
J=lambda p:json.loads(p.read_text());pin=lambda b:dict(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
canonical=lambda v:pin(json.dumps(v,sort_keys=True,separators=(',',':')).encode())
assert not (P/'provider.patch').exists()
old=J(V/'source-contract.json');before={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
assert canonical(before)==old['full_candidate_source_identity']
cpp=J(W/'cpp-build-result.json');assert cpp['passed'] and cpp['new_CPP_compilations']==1 and cpp['new_relocatable_links']==1
probe=J(W/'cpp/run/stdout.log');assert probe['passed'] and probe['compared_words']==44452800 and probe['tag_words']==8890560
assert J(W/'cpp/run/owner.json')['passed'] and J(W/'cpp-code-inspection.json')['passed']
objects=old['native_objects'][:7]+[cpp['merged_adapter']];assert len(objects)==8
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
texts={};name='src/gallium/targets/lavapipe/meson.build';text=(S/name).read_text()
a,b=old['native_objects'][-1],objects[-1]
assert text.count(str(R/a['path']))==1;texts[name]=text.replace(str(R/a['path']),str(R/b['path']),1)
for name,a,b in [('src/gallium/drivers/llvmpipe/lp_screen.c','Sirius-native-cpu-normalization-tags-v34','Sirius-native-cpu-predicated-refinement-v35'),('src/gallium/frontends/lavapipe/lvp_device.c','SIR-CPU-25.2.8AJ','SIR-CPU-25.2.8AK')]:
 text=(S/name).read_text();assert text.count(a)==1;texts[name]=text.replace(a,b,1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
 original=(S/name).read_text();patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,fromfile='a/'+name,tofile='b/'+name))
 after[name]=pin(text.encode());modified[name]=dict(basis=before[name],candidate=after[name])
inputs=dict(old['source_inputs'])
inputs.update(J(W/'prepared-cpp.json')['inputs'])
for row in objects:inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
for p in [V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',W/'cpp-build-result.json',W/'cpp-code-inspection.json',W/'cpp/run/owner.json',W/'cpp/run/stdout.log',W/'cpp-source-proof.json',Path(__file__).resolve()]:inputs[str(p.relative_to(R))]=pin(p.read_bytes())
for name,v in inputs.items():assert pin((R/name).read_bytes())==v,name
contract=dict(basis_source_identity=old['full_candidate_source_identity'],full_candidate_source_identity=canonical(after),modified_files=modified,patch=pin(patch.encode()),
 native_objects=objects,source_inputs=inputs,new_native_compilations=1,new_ISPC_compilations=0,new_CPP_compilations=1,reused_raw_objects=3,
 guard_shared_byte_authority='out/native-cpu-active-prefix/native_packet.hpp',
 scope='AK restores six empty-work returns in original C++ Divide/Sqrt only. Exact source inverse restores AJ. Original eligible identity/refusal/invalid-payload/zero selection is returned when no refinement lane remains; every nonempty work mask retains the original seed, three ordered corrections, final residual, error, pack and refusal. Same complete five words, tag predicate/write protocol, all raw3/scalar3/pure admission/original adapter,9380tuples/SPVs/layouts, eligibility, FP/ISA, counts/cache/governors, fields/RHS/trials/quartics/order/rollback/scheduler. One new helper compilation, existing fixture driver and seven provider objects reused. Fresh whole-bank/tag fixture passes; actual provider joins/full38/C8/affected whole-path checks still required. No gain/science/release claim.')
P.mkdir(exist_ok=True);(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps(dict(modified_provider_files=3,new_helper_compilations=1,reused_raw_objects=3,shared_Mesa_modified=False)))
