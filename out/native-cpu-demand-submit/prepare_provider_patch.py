"""Keep existing dependency-driven queue submission, reusing all native objects."""
from pathlib import Path
import difflib
import hashlib
import json
import subprocess

R=Path.cwd();W=Path(__file__).resolve().parent;P=W/'provider'
S=R/'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
V=R/'out/native-cpu-masked-arithmetic/provider'
doc=lambda p:json.loads(p.read_text())
pin=lambda data:{'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
canonical=lambda d:pin(json.dumps(d,sort_keys=True,separators=(',',':')).encode())
assert not (P/'provider.patch').exists()
failed=R/'out/native-cpu-masked-arithmetic-cohort/cold-frame'
owner=doc(failed/'frame/owner.json');assert owner['status']=='terminal' and not owner['passed']
assert owner['stop_reason']=='diagnostic_timeout' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
manifest=R/'attestations/qualification'/owner['source']['revision'][:8]/'native-cpu-masked-arithmetic-frame-timeout/MANIFEST.json'
assert not doc(manifest)['accepted'] and not doc(failed/'result.json')['accepted']
assert subprocess.check_output(['git','show','@{upstream}:'+str(manifest.relative_to(R))])==manifest.read_bytes()
old=doc(V/'source-contract.json')
before={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
assert canonical(before)==old['full_candidate_source_identity']
objects=old['native_objects'];assert len(objects)==8
for row in objects:assert pin((R/row['path']).read_bytes())=={k:row[k] for k in ('bytes','sha256')}
texts={}
name='src/gallium/frontends/lavapipe/lvp_device.c';original=(S/name).read_text()
force='''   result = vk_queue_enable_submit_thread(&queue->vk);
   if (result != VK_SUCCESS) {
      vk_queue_finish(&queue->vk);
      return result;
   }
'''
replacement='''#if !LP_VULKAN_IEEE32_ENABLED
'''+force+'''#endif
   /* IEEE32 keeps Mesa's on-demand submit mode. Ready submissions execute
    * through the same driver callback; unresolved dependencies start the
    * existing submit thread in vk_queue_submit(). Queue order, sync payloads,
    * driver FP save/restore, execution and completion remain unchanged.
    */
'''
assert original.count(force)==1 and original.count('SIR-CPU-25.2.8AE')==1
texts[name]=original.replace(force,replacement,1).replace('SIR-CPU-25.2.8AE','SIR-CPU-25.2.8AF',1)
assert texts[name].replace(replacement,force,1).replace('SIR-CPU-25.2.8AF','SIR-CPU-25.2.8AE',1)==original
name='src/gallium/drivers/llvmpipe/lp_screen.c';text=(S/name).read_text()
assert text.count('Sirius-native-cpu-masked-arithmetic-v29')==1
texts[name]=text.replace('Sirius-native-cpu-masked-arithmetic-v29','Sirius-native-cpu-demand-submit-v30',1)
runtime={name:(S/name).read_text() for name in ('src/vulkan/runtime/vk_device.c','src/vulkan/runtime/vk_queue.c')}
assert 'device->submit_mode = VK_QUEUE_SUBMIT_MODE_THREADED_ON_DEMAND;' in runtime['src/vulkan/runtime/vk_device.c']
q=runtime['src/vulkan/runtime/vk_queue.c']
for fragment in ('VK_SYNC_WAIT_PENDING, 0','result = vk_queue_enable_submit_thread(queue);','result = vk_queue_submit_final(queue, submit);','vk_queue_submit_move_binary_waits_to_temps(device, submit, false)'):
 assert fragment in q
patch='';after=dict(before);modified={}
for name,text in texts.items():
 original=(S/name).read_text()
 patch+=''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),n=3,fromfile='a/'+name,tofile='b/'+name))
 after[name]=pin(text.encode());modified[name]={'basis':before[name],'candidate':after[name]}
inputs=dict(old['source_inputs'])
for p in [V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',failed/'frame/owner.json',failed/'result.json',manifest,W/'prepare_provider_patch.py',*[S/name for name in runtime]]:
 inputs[str(p.relative_to(R))]=pin(p.read_bytes())
for name,expected in inputs.items():assert pin((R/name).read_bytes())==expected,name
contract={'basis_source_identity':old['full_candidate_source_identity'],'full_candidate_source_identity':canonical(after),
 'modified_files':modified,'patch':pin(patch.encode()),'native_objects':objects,'source_inputs':inputs,
 'new_native_compilations':0,'guard_shared_byte_authority':'out/native-cpu-active-prefix/native_packet.hpp',
 'scope':'AF removes unconditional queue-thread startup only in the private IEEE32 build, retaining Mesa existing on-demand mode and all original unresolved-dependency thread fallback, queue ordering, binary/timeline sync, failure, execution and completion paths. Same driver callback establishes/restores FP state regardless of caller. Non-IEEE builds retain original forced startup. Reuse every exact AE native object, code and fresh full38 evidence; no native math, ISA, eligibility, tuples, counts, shader, provider options or budgets change. This source argument is not machine-checked proof. Fresh provider code joins/C8/original frame required; broader Z evidence remains predecessor only. No speed, scientific or release claim.'}
P.mkdir(exist_ok=True);(P/'provider.patch').write_text(patch);(P/'source-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
assert before=={str(p.relative_to(S)):pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_files':2,'new_native_compilations':0,'shared_Mesa_modified':False}))
