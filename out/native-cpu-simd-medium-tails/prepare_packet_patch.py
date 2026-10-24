"""Emit an unapplied packet-policy patch and a closed row-partition proof."""
from pathlib import Path
import difflib
import hashlib
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
mesa = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
driver = 'src/gallium/drivers/llvmpipe/'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        path.write_bytes(data)


def replace_once(text, old, new):
    assert text.count(old) == 1, old[:120]
    return text.replace(old, new, 1)


base_hashes = {'lp_state_cs.c': '0f242787a0d7ba08ca35157fccfbc0bde12c6538f8101669a5c59a34e7f166c7',
               'lp_native_cpu_bridge.h': 'ffa014c8d06e66c33935477fa11a13b775b9cc479860a3be5d6943adc7c3568f'}
basis = {}
inputs = {}
for name, expected in base_hashes.items():
    path = mesa / driver / name
    data = path.read_bytes()
    assert pin(data)['sha256'] == expected, name
    basis[name] = data.decode()
    inputs[str(path.relative_to(root))] = pin(data)

changes = {
    'lp_state_cs.c': [
        ('''      /* Full packets cover eight rows; each remainder row owns one task. */
      unsigned full_packets = job_info->grid_size[0] / 8;
      grid_z = grid_y = 0;
      grid_x = (unsigned)iter_idx < full_packets ? (unsigned)iter_idx * 8 :
         full_packets * 8 + ((unsigned)iter_idx - full_packets);
''', '''      /* Full packets cover eight rows; medium tails share one masked task.
       * One-to-three-row tails retain independent original scalar tasks. */
      unsigned full_packets = job_info->grid_size[0] / 8;
      unsigned remainder = job_info->grid_size[0] % 8;
      grid_z = grid_y = 0;
      grid_x = (unsigned)iter_idx < full_packets ? (unsigned)iter_idx * 8 :
         full_packets * 8 + (remainder >= 4 ? 0 : (unsigned)iter_idx - full_packets);
'''),
        ('''      scheduled_tasks = job_info.grid_size[0] / 8 + job_info.grid_size[0] % 8;
''', '''      unsigned remainder = job_info.grid_size[0] % 8;
      scheduled_tasks = job_info.grid_size[0] / 8 + (remainder >= 4 ? 1 : remainder);
'''),
    ],
    'lp_native_cpu_bridge.h': [
        ('''      groups = x < (job->grid_size[0] / 8) * 8 ? 8 : 1;
''', '''      groups = job->grid_size[0] - x >= 4 ? MIN2(8, job->grid_size[0] - x) : 1;
'''),
    ],
}
new_sources = {}
patch, modified = '', {}
for name, replacements in changes.items():
    new = basis[name]
    for old, replacement in replacements:
        new = replace_once(new, old, replacement)
    restored = new
    for old, replacement in reversed(replacements):
        restored = replace_once(restored, replacement, old)
    assert restored == basis[name]
    relative = driver + name
    new_sources[name] = new
    modified[relative] = {'basis': pin(basis[name].encode()), 'candidate': pin(new.encode())}
    patch += ''.join(difflib.unified_diff(basis[name].splitlines(True), new.splitlines(True),
                                         n=3, fromfile='a/' + relative, tofile='b/' + relative))

i = root / 'out/native-cpu-simd-packet-policy'
old_adapter = (i / 'rows_adapter.cpp').read_bytes()
assert pin(old_adapter)['sha256'] == '57747ee37daa040d25b1e6f9411affb07fb399ba0003c921aca313810a10a581'
old = b'if (groups != 8U || !sirius_native_simd_supported())'
new = b'if (groups < 4U || groups > 8U || !sirius_native_simd_supported())'
adapter = (work / 'rows_adapter.cpp').read_bytes()
assert old_adapter.count(old) == adapter.count(new) == 3
assert adapter.replace(new, old) == old_adapter
contract_path = root / 'out/native-cpu-nir-dead-progress/provider/source-contract.json'
original = json.loads(contract_path.read_text())
objects = original['native_objects']
assert len(objects) == 8
derivation_path = root / 'out/native-cpu-simd-endpoint-regions/ispc-derivation.json'
derivation = json.loads(derivation_path.read_text())
assert derivation['instructions'] == 9380
for slug in ('transport', 'endpoint', 'dopri_phase'):
    item = derivation['stages'][slug]
    path = root / item['derived_path']
    data = path.read_bytes()
    assert pin(data) == item['derived']
    text = data.decode()
    # The already compiled raw function masks both execution and publication;
    # no source or arithmetic changes are needed to admit four-to-seven rows.
    assert 'if (programIndex < groups) ComputeOwnedRow_0(id, invocation, &context);' in text
    assert '    if (programIndex < groups) {\n        for (uniform uint32_t column' in text
    assert item['logical_invocations_per_active_row'] == 64
    inputs[str(path.relative_to(root))] = pin(data)
for path in (i / 'rows_adapter.cpp', i / 'adapter-compile-result.json', contract_path, derivation_path,
             root / 'out/native-cpu-active-prefix/native_packet.hpp', work / 'rows_adapter.cpp', work / 'rows_bridge.h'):
    inputs[str(path.relative_to(root))] = pin(path.read_bytes())

# Euclidean division is a closed proof for every N, without enumerating a
# workload or sampling shapes: N=8q+r, 0<=r<8. The eight residue cases below
# partition the entire tail interval after q disjoint eight-row intervals.
residues = []
for remainder in range(8):
    widths = [remainder] if remainder >= 4 else [1] * remainder
    assert sum(widths) == remainder
    offset, intervals = 0, []
    for width in widths:
        assert width == 1 or 4 <= width <= 7
        intervals.append({'start': '8q+' + str(offset), 'width': width})
        offset += width
    assert offset == remainder
    residues.append({'r': remainder, 'tail_tasks': len(widths), 'tail_intervals': intervals,
                     'total_tasks': 'q+' + str(len(widths))})
proof = {
    'patch': pin(patch.encode()), 'changed_files': modified, 'source_inputs': inputs,
    'adapter_only_three_group_predicates_changed': True,
    'seven_reused_native_objects': objects[:7], 'new_native_objects_required': 1,
    'raw_ISPC_sources_and_all9380calls_unchanged': True,
    'partition': {
        'domain': 'Every integer N with1<=N<=65536; immutable B admission additionally proves N<=capacity<=65536.',
        'division': 'Unique q=floor(N/8),r=N%8 withN=8q+r and0<=r<8.',
        'full_tasks': 'For0<=i<q, x=8i,groups8, interval[8i,8i+8). Their disjoint union is[0,8q).',
        'medium_tail': 'For4<=r<=7, exactly one taski=q hasx=8q,groupsr, interval[8q,N).',
        'small_tail': 'For0<=r<=3, taski=q+j,0<=j<r, hasx=8q+j,groups1. These adjacent singletons partition[8q,N).',
        'coverage': 'The full and tail intervals are adjacent, pairwise disjoint and have total length8q+r=N. Every original row executes exactly once and no row at or aboveN executes.',
        'C_integer_safety': 'q<=8192 and0<=r<=7; taskcount<=8194 forN<=65536 (largest singleton tail atN=65531). All iterator arithmetic fits unsigned/signed task ABI. Existing coarse bounds remain unchanged.',
        'residues': residues,
        'CPU_false': 'native_simd_packets=false retains the original N scalar tasksx=i,groups1 for0<=i<N. The adapter independently forwards every call to unchanged D scalar when the existing CPU/v3 guard is false.',
        'query_count': 'Original num_tasks=N andblock64 are untouched; cs_invocations remainsuint64_t(N)*64, independently of scheduled_tasks.',
        'zero_refusal_fallback': 'Existing separate zero-header classification, positive B/F pure admission and ordinary LLVM key/task/fallback branches are unchanged.',
    },
    'inactive_lane_proof': 'Raw L exports already mask ComputeOwnedRow and final output copy withprogramIndex<groups. Private context/staged output remain per-row; no unexecuted lane reads an input row or publishes output. B first_x+groups<=activegrid<=capacity protects all active rows.',
    'limits': 'Source-only medium-tail utilization correction motivated by saved frame cardinality histograms. The existing samples do not establish the full-frame cause or a speedup. No shader/compiler/provider/workload execution or performance/scientific/release claim.',
    'patch_applied': False,
}
publish(work / 'packet-policy.patch', patch.encode())
publish(work / 'partition-proof.json', (json.dumps(proof, indent=2) + '\n').encode())
assert all((mesa / driver / name).read_text() == basis[name] for name in basis)
print(json.dumps({'source_patch_written': True, 'all_N_partition_proved': True,
                  'seven_native_objects_reused': 7, 'new_adapter_required': 1,
                  'patch_applied': False, 'compiled_or_executed': False}))
