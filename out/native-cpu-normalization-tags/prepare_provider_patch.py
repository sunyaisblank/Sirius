"""Replace native raw3/adapter and private identities, retaining exact Mesa source."""
from pathlib import Path
import difflib,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider'
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8';V=R/'out/native-cpu-direct-compound/provider'
J=lambda p:json.loads(p.read_text());pin=lambda b:{'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
canonical=lambda v:pin(json.dumps(v,sort_keys=True,separators=(',',':')).encode())
assert not (P/'provider.patch').exists()
old=J(V/'source-contract.json');before={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()};assert canonical(before)==old['full_candidate_source_identity']
compiled=J(W/'ispc-dense-compile-result.json');assert compiled['all_objects_compiled'] and compiled['new_compilations']==4 and compiled['new_raw_compilations']==3 and compiled['new_predicate_compilations']==1
assert J(W/'dense-code-inspection.json')['passed'] and J(W/'object-code-inspection-dense.json')['objects']
cpp=J(W/'cpp-build-result.json');assert cpp['passed'] and cpp['new_relocatable_links']==1
probe=J(W/'cpp/run/stdout.log');assert probe['passed'] and probe['compared_words']==44452800 and probe['tag_words']==8890560 and J(W/'cpp/run/owner.json')['passed']
assert J(W/'predicate/result.json')['passed'] and J(W/'predicate/run/owner.json')['passed']
objects=old['native_objects'][:4]+[row['object'] for row in compiled['objects']]+[cpp['merged_adapter']];assert len(objects)==8
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
texts={};name='src/gallium/targets/lavapipe/meson.build';text=(S/name).read_text()
for a,b in zip(old['native_objects'][4:],objects[4:]):
 assert text.count(str(R/a['path']))==1;text=text.replace(str(R/a['path']),str(R/b['path']),1)
texts[name]=text
for name,a,b in (('src/gallium/drivers/llvmpipe/lp_screen.c','Sirius-native-cpu-direct-compound-v33','Sirius-native-cpu-normalization-tags-v34'),('src/gallium/frontends/lavapipe/lvp_device.c','SIR-CPU-25.2.8AI','SIR-CPU-25.2.8AJ')):
 text=(S/name).read_text();assert text.count(a)==1;texts[name]=text.replace(a,b,1)
patch='';after=dict(before);modified={}
for name,text in texts.items():
 original=(S/name).read_text();patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,fromfile='a/'+name,tofile='b/'+name));after[name]=pin(text.encode());modified[name]={'basis':before[name],'candidate':after[name]}
# Current compiler/source/object dependencies are sealed directly. Historical
# unchanged native4 provenance is bound by its source contract and linked proof,
# rather than recursively importing every expired predecessor's diagnostic file.
inputs=dict(compiled['source_inputs'])
for row in objects:inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
for p in (V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',W/'cpp-build-result.json',W/'cpp-code-inspection.json',W/'cpp/run/owner.json',W/'cpp/run/stdout.log',W/'ispc-dense-compile-result.json',W/'dense-code-inspection.json',W/'object-code-inspection-dense.json',W/'compile/ispc-dense/owner.json',W/'predicate/result.json',W/'predicate/link-result.json',W/'predicate/run/owner.json',W/'predicate/link/owner.json',W/'normalization_tag.isph',W/'inspect_tagged_objects.py',W/'inspect_tagged_dense.py',W/'object-inspection-correction.json',Path(__file__).resolve()):inputs[str(p.relative_to(R))]=pin(p.read_bytes())
for name,v in inputs.items():assert pin((R/name).read_bytes())==v,name
contract={'basis_source_identity':old['full_candidate_source_identity'],'full_candidate_source_identity':canonical(after),'modified_files':modified,'patch':pin(patch.encode()),'native_objects':objects,'source_inputs':inputs,'new_native_compilations':6,'new_ISPC_compilations':4,'new_CPP_compilations':2,'reused_raw_objects':0,'guard_shared_byte_authority':'out/native-cpu-active-prefix/native_packet.hpp',
 'scope':'AJ caches only the unchanged normalized-input predicate alongside each complete five-word private register. The complete context initializes metadata; both generic/uniform partial writes invalidate it and complete writes/direct CPP stores recompute it from the same actual words. Tagged Add/Subtract/Multiply consume only current exact tags for valid lanes, with the original high2^54 Multiply guard and complete fallback. Original AI CPP source is an exact prefix; Divide/Sqrt and all formulas/operation order/normalization/error/refusal remain. No presumed input/postcondition normalization, persistent/input/shader-cache change, extra published limb, stronger eligibility, altered FP/ISA/precision/requests/fields/stages/trials/quartics/counts/row order/rollback/scheduler/queue/governors. All9380literal tuples and public layouts remain. Fresh370440call/44452800payload-word/8890560tag-word original-scalar/std::nextafter comparison plus exact shared ISPC predicate fixture passed; actual linked targets/full38/C8/affected control5/sole original frame required. H host/frame reused unchanged. No performance/science/release acceptance.'}
P.mkdir(exist_ok=True);(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_provider_files':3,'new_raw_objects':3,'shared_Mesa_modified':False}))
