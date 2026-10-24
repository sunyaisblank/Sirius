"""Derive the Transport worker-occupancy policy without changing Mesa."""
from pathlib import Path
import difflib
import hashlib
import json
import re

root = Path.cwd()
work = Path(__file__).resolve().parent
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
name = 'src/gallium/drivers/llvmpipe/lp_state_cs.c'
source = (source_root / name).read_text()


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        path.write_bytes(data)


anchor = '#include "lp_native_cpu_bridge.h"\n'
helper = '''
/* Pure policy only: retain serial-pool packets; expose one row per worker
 * when all Transport rows fit the existing pool. No pool setting changes.
 */
static bool ATTRIBUTE_NOINLINE
sirius_native_transport_simd_policy(unsigned rows, unsigned workers)
{
   return workers == 0 || rows > workers;
}
'''
old = '         job_info.native_simd_packets = !native_zero_active && sirius_native_simd_supported();'
replacement = '''         job_info.native_simd_packets = !native_zero_active && sirius_native_simd_supported() &&
            (job_info.native_cpu_stage != 1 ||
             sirius_native_transport_simd_policy(job_info.grid_size[0], screen->cs_tpool->num_threads));'''
assert source.count(anchor) == source.count(old) == 1
assert '#include "lp_cs_tpool.h"' in source
new = source.replace(anchor, anchor + helper, 1).replace(old, replacement, 1)
assert new.replace(replacement, old, 1).replace(anchor + helper, anchor, 1) == source
pool_name = 'src/gallium/drivers/llvmpipe/lp_cs_tpool.c'
pool = (source_root / pool_name).read_text()
assert pool.count('pool->num_threads = num_threads;') == 1
assert len(re.findall(r'pool->num_threads\s*=(?!=)', pool)) == 1
assert 'if (pool->num_threads == 0)' in pool
assert 'task->iter_per_thread = num_iters / pool->num_threads;' in pool
assert 'task->iter_remainder = num_iters % pool->num_threads;' in pool
bridge_name = 'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h'
bridge = (source_root / bridge_name).read_text()
assert 'case 1:\n      return sirius_native_transport_rows_try_packet' in bridge
assert 'if (job->native_simd_packets)\n      groups = x < (job->grid_size[0] / 8) * 8 ? 8 : 1;' in bridge
assert 'if (job_info->native_simd_packets)' in source
assert 'scheduled_tasks = job_info.grid_size[0] / 8 + job_info.grid_size[0] % 8;' in source
assert '(uint64_t)num_tasks * info->block[0] * info->block[1] * info->block[2]' in source
new_text = {name: new}
for filename, before, after in (
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8M', 'SIR-CPU-25.2.8Q'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-nir-dead-progress-v13',
     'Sirius-native-cpu-transport-worker-occupancy-v15'),
):
    text = (source_root / filename).read_text()
    assert text.count(before) == 1
    new_text[filename] = text.replace(before, after, 1)
patch = ''
modified = {}
for filename, text in new_text.items():
    before = (source_root / filename).read_text()
    patch += ''.join(difflib.unified_diff(before.splitlines(True), text.splitlines(True),
        fromfile='a/' + filename, tofile='b/' + filename, n=3))
    modified[filename] = {'basis': pin(before.encode()), 'candidate': pin(text.encode())}
proof = {
    'changed_files': modified, 'patch': pin(patch.encode()), 'patch_applied': False,
    'policy_helper': 'sirius_native_transport_simd_policy',
    'policy': 'positive && CPU_ISA && (stage != Transport || workers == 0 || grid_x > workers)',
    'policy_source_inverse_exact': True, 'new_native_compilations': 0,
    'source_inputs': {str((source_root / filename).relative_to(root)):
                     pin((source_root / filename).read_bytes())
                     for filename in (*new_text, pool_name, 'src/gallium/drivers/llvmpipe/lp_cs_tpool.h', bridge_name)},
    'partition_proof': {
        'domain': 'All admitted original grids 1 <= N <= 65536 and every unsigned pool worker count W.',
        'Transport_CPU_true_W_positive_N_le_W': 'The existing nonpacket native branch queues N tasks, IDs 0..N-1, groups1. Each original row runs once through the unchanged F scalar execution path.',
        'other_positive_jobs': 'E/DP, Transport W0, Transport N>W retain M exactly: q=floor(N/8), r=N%8; q+r tasks map 8*i for i<q (groups8), then 8*q+i-q (groups1). Disjoint intervals partition [0,N).',
        'CPU_false': 'Original F scalar N-task path for all stages, unchanged.',
        'refusal_zero_active_and_other_modules': 'Policy not evaluated. All admission, zero-active no-op, ordinary intergroup key/JIT, resource/dirty ordering and failure paths remain M.',
        'queries': 'num_tasks remains original N; both paths add N*64 logical invocations. Worker FP scope and original guards remain unchanged.',
        'pool_lifetime': 'num_threads assigned once after creation/startup failures, never changed after publication; no environment/worker/governor setting altered.'},
    'cost_scope': 'M Transport12 is one8-row SIMD task plus four scalar tasks; Q can expose twelve one-row tasks when the existing pool has at least twelve workers. No proven latency crossover, frame bottleneck attribution or speedup. This is a topology-based implementation policy, not a changed scientific gate or threshold search.',
    'gate_scope': 'All original C8 32 readbacks/51 markers/54 joins and original ABBA22/Decimal thresholds remain unchanged. Only extra passive implementation-route expectations adapt; historical Transport SIMD execution is not fresh Q execution evidence.'}
publish(work / 'provider.patch', patch.encode())
publish(work / 'source-proof.json', (json.dumps(proof, indent=2) + '\n').encode())
assert (source_root / name).read_text() == source
print(json.dumps({'unapplied_patch': True, 'changed_files': 3, 'native_compilations': 0,
                  'logical_rows_and_queries_unchanged': True}))
