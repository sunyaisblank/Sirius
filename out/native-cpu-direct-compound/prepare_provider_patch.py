"""Replace native raw3/adapter and private identities, retaining exact Mesa source."""
from pathlib import Path
import difflib,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider'
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8';V=R/'out/native-cpu-normalized-arithmetic/provider'
J=lambda p:json.loads(p.read_text());pin=lambda b:{'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
canonical=lambda v:pin(json.dumps(v,sort_keys=True,separators=(',',':')).encode())
assert not (P/'provider.patch').exists()
old=J(V/'source-contract.json');before={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()};assert canonical(before)==old['full_candidate_source_identity']
compiled=J(W/'ispc-dense-compile-result.json');assert compiled['all_objects_compiled'] and compiled['new_compilations']==3
assert J(W/'dense-code-inspection.json')['passed'] and J(W/'object-code-inspection-dense.json')['objects']
cpp=J(W/'cpp-build-result.json');assert cpp['passed'] and cpp['new_relocatable_links']==1
probe=J(W/'cpp/run/stdout.log');assert probe['passed'] and probe['compared_words']==22226400 and J(W/'cpp/run/owner.json')['passed']
objects=old['native_objects'][:4]+[row['object'] for row in compiled['objects']]+[cpp['merged_adapter']];assert len(objects)==8
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
texts={};name='src/gallium/targets/lavapipe/meson.build';text=(S/name).read_text()
for a,b in zip(old['native_objects'][4:],objects[4:]):
 assert text.count(str(R/a['path']))==1;text=text.replace(str(R/a['path']),str(R/b['path']),1)
texts[name]=text
for name,a,b in (('src/gallium/drivers/llvmpipe/lp_screen.c','Sirius-native-cpu-normalized-arithmetic-v32','Sirius-native-cpu-direct-compound-v33'),('src/gallium/frontends/lavapipe/lvp_device.c','SIR-CPU-25.2.8AH','SIR-CPU-25.2.8AI')):
 text=(S/name).read_text();assert text.count(a)==1;texts[name]=text.replace(a,b,1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
 original=(S/name).read_text();patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,fromfile='a/'+name,tofile='b/'+name));after[name]=pin(text.encode());modified[name]={'basis':before[name],'candidate':after[name]}
# Current compiler/source/object dependencies are sealed directly. Historical
# unchanged native4 provenance is bound by its source contract and linked proof,
# rather than recursively importing every expired predecessor's diagnostic file.
inputs=dict(compiled['source_inputs'])
for row in objects:inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
for p in (V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',W/'cpp-build-result.json',W/'cpp-code-inspection.json',W/'cpp/run/owner.json',W/'cpp/run/stdout.log',W/'ispc-dense-compile-result.json',W/'dense-code-inspection.json',W/'object-code-inspection-dense.json',W/'compile/ispc-dense/owner.json',Path(__file__).resolve()):inputs[str(p.relative_to(R))]=pin(p.read_bytes())
for name,v in inputs.items():assert pin((R/name).read_bytes())==v,name
contract={'basis_source_identity':old['full_candidate_source_identity'],'full_candidate_source_identity':canonical(after),'modified_files':modified,'patch':pin(patch.encode()),'native_objects':objects,'source_inputs':inputs,'new_native_compilations':5,'new_ISPC_compilations':3,'new_CPP_compilations':2,'reused_raw_objects':0,'guard_shared_byte_authority':'out/native-cpu-active-prefix/native_packet.hpp',
 'scope':'AI changes static opcode3/5/6 dispatch to the same direct private bank already used for Add/Multiply, removing zeroed temporary banks, operand/result repacking and returned tuples at936sites. All9380ordered tuples and original generic/value helpers remain, with whole-source inverse proof. New subtract entry is exactly unchanged AH Add(a,Negate(b)); all existing AH source is an exact prefix. All arithmetic/final checks/FP/ISA/eligibility/requests/fields/trials/quartics/row order/refusal/rollback/scheduler/queue/governors remain. Five-operation49case22226400word complete-bank comparison passed; actual linked targets/full38/C8/affected original control5 and sole original frame still required. H host/frame are reused unchanged. Source/finite compatibility only, no performance/science/release acceptance.'}
P.mkdir(exist_ok=True);(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_provider_files':3,'new_raw_objects':3,'shared_Mesa_modified':False}))
