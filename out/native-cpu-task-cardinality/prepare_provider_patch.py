"""Use exactly available work to notify workers; directly run one native task."""
from pathlib import Path
import difflib
import hashlib
import json
import re

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
X = R / 'out/native-cpu-predicated-arithmetic/provider'
doc = lambda p: json.loads(p.read_text())
pin = lambda b: {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
canonical = lambda value: pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())

def worker(text):
    match = re.search(r'static void\ncs_exec_fn\([^\n]+\)\n\{', text); assert match
    end = match.end(); depth = 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}'); end += 1
    return text[match.start():end]

assert not (P / 'provider.patch').exists()
old = doc(X / 'source-contract.json')
before = {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
assert canonical(before) == old['full_candidate_source_identity']
texts = {}
name = 'src/gallium/drivers/llvmpipe/lp_cs_tpool.c'; text = (S / name).read_text()
a = '''   list_addtail(&task->list, &pool->workqueue);

   cnd_broadcast(&pool->new_work);
   mtx_unlock(&pool->m);'''
b = '''   list_addtail(&task->list, &pool->workqueue);

   /* A short task has at most one iteration for each notified worker.
    * Already-running workers also recheck this queue under the same mutex.
    * Keep broadcast when every pool worker can have work; otherwise avoid
    * waking workers for which this task has no iteration.
    */
   if ((unsigned)num_iters >= pool->num_threads) {
      cnd_broadcast(&pool->new_work);
   } else {
      for (unsigned i = 0; i < (unsigned)num_iters; i++)
         cnd_signal(&pool->new_work);
   }
   mtx_unlock(&pool->m);'''
assert text.count(a) == 1; texts[name] = text.replace(a, b, 1)
assert texts[name].count('cnd_broadcast(&pool->new_work)') == 2  # Work and shutdown.
name = 'src/gallium/drivers/llvmpipe/lp_state_cs.c'; text = (S / name).read_text()
a = '''   if (scheduled_tasks) {
      struct lp_cs_tpool_task *task;'''
b = '''   if (job_info.native_preflight_passed && scheduled_tasks == 1) {
      /* This dispatch is synchronous and native admission proves no shared
       * allocation. Use the same worker/IEEE32 scope directly for its sole
       * task, without a queue wake and immediate wait. Ordinary work retains
       * its existing pool route; logical invocation accounting stays below.
       */
      struct lp_cs_local_mem native_lmem = {0};
      cs_exec_fn(&job_info, 0, &native_lmem);
      FREE(native_lmem.local_mem_ptr);
   } else if (scheduled_tasks) {
      struct lp_cs_tpool_task *task;'''
assert text.count(a) == 1; candidate = text.replace(a, b, 1)
assert worker(text) == worker(candidate)
texts[name] = candidate
for name, a, b in (
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-predicated-arithmetic-v22', 'Sirius-native-cpu-task-cardinality-v23'),
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8X', 'SIR-CPU-25.2.8Y')):
    text = (S / name).read_text(); assert text.count(a) == 1; texts[name] = text.replace(a, b, 1)
patch = ''; after = dict(before); modified = {}
for name, text in texts.items():
    original = (S / name).read_text()
    patch += ''.join(difflib.unified_diff(original.splitlines(True), text.splitlines(True), n=3, fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(text.encode()); modified[name] = {'basis': before[name], 'candidate': after[name]}
objects = old['native_objects']; assert len(objects) == 8
for row in objects:
    assert pin((R / row['path']).read_bytes()) == {k: row[k] for k in ('bytes', 'sha256')}
inputs = dict(old['source_inputs'])
fp_path = str((S / 'src/gallium/drivers/llvmpipe/lp_state_cs.c').relative_to(R))
assert inputs[fp_path] == before['src/gallium/drivers/llvmpipe/lp_state_cs.c']
for path in (X / 'source-contract.json', X / 'provider-build-result.json', X / 'linked-native-proof.json',
             R / 'out/native-cpu-predicated-arithmetic-compatibility/result.json', W / 'prepare_provider_patch.py'):
    inputs[str(path.relative_to(R))] = pin(path.read_bytes())
for name, expected in inputs.items():
    assert pin((R / name).read_bytes()) == expected, name
inputs[fp_path] = after['src/gallium/drivers/llvmpipe/lp_state_cs.c']
contract = {'basis_source_identity': old['full_candidate_source_identity'],
    'full_candidate_source_identity': canonical(after), 'modified_files': modified, 'patch': pin(patch.encode()),
    'native_objects': objects, 'source_inputs': inputs, 'new_native_compilations': 0,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'FP_ownership_source_transition': {'path': fp_path, 'basis': before['src/gallium/drivers/llvmpipe/lp_state_cs.c'],
        'candidate': inputs[fp_path], 'identical_cs_exec_fn_body': pin(worker(texts['src/gallium/drivers/llvmpipe/lp_state_cs.c']).encode()),
        'reason': 'X derivation pinned this file only to establish runtime FP ownership. The exact worker/FP body is unchanged; only launch scheduling changes. No arithmetic compiler input or native object changes.'},
    'scope': 'Y keeps all8 native object bytes, complete programs/SPIRV, packet/row choices, pool thread count, numerical/FP/eligibility/field/count/query/fallback/rollback contracts. Single positively admitted native task uses the identical worker synchronously outside cs_mutex; short task notifications follow actual cardinality, with existing queue mutex/claim/finish protocol and shutdown broadcast retained. No native/host/fixture rebuild, cache-policy/governor/threshold tuning or speed/qualification claim.',
    'reasoning': [
        'A queued positive task is visible before notification under pool->m. Each notified worker claims unique iterations under that mutex. Existing running workers recheck the same queue before waiting; when nobody waits, the visible queue still supplies work. Larger tasks retain broadcast. Partition/finish counters and wait predicates are unchanged.',
        'The native direct branch requires completed pure admission and exactly one original scheduled task. It executes the same cs_exec_fn, same index0 and zeroed local-memory state, before identical logical query accounting. Native shape admission excludes shared/payload allocations. The existing zero-worker route already establishes synchronous worker execution semantics; ordinary and zero-active paths remain.',
        'No assumption that the single cumulative OS snapshot attributes thread CPU time to native math or proves a causal share. The structural premise is the unconditional source broadcast and queue/wait for one admitted task. Original C8/ABBA22+Decimal/all8/frame remain decisive.'
    ]}
P.mkdir(exist_ok=True)
(P / 'provider.patch').write_text(patch)
(P / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_provider_files': 4, 'reused_native_objects': 8, 'new_native_compilations': 0, 'shared_Mesa_modified': False}))
