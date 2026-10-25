"""Replace row SIMD with within-row field SIMD, retaining logical scheduling."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'


def doc(path):
    return json.loads(path.read_text())


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(value):
    return pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


assert not (P / 'provider.patch').exists()
old = doc(R / 'out/native-cpu-zero-term/provider/source-contract.json')
before = {str(path.relative_to(S)): pin(path.read_bytes()) for path in S.rglob('*') if path.is_file()}
assert canonical(before) == old['full_candidate_source_identity']
objects = old['native_objects'][:4]
inputs = {}
for kind in ('cpp', 'ispc', 'adapter'):
    name = kind + '-compile-result.json'
    compiled = doc(W / name); owner = doc(W / 'compile' / kind / 'owner.json')
    assert owner['passed'] and owner['whole_inputs_unchanged'] and owner['source_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['exit'] == 0 for row in compiled['objects'])
    if kind != 'adapter':
        assert compiled['all_objects_compiled'] and compiled['new_compilations'] == 3
    inputs.update(compiled['source_inputs'])
    objects += [row['object'] for row in compiled['objects']]
    for path in (W / name, W / 'compile' / kind / 'owner.json'):
        inputs[str(path.relative_to(R))] = pin(path.read_bytes())
assert len(objects) == 11
isa = doc(W / 'object-code-inspection.json')
for slug in ('transport', 'endpoint', 'dopri_phase'):
    row = isa['objects'][slug]
    assert not row['forbidden_fma_fp_control_x87_found'] and not row['forbidden_non_v3_extensions'] and not row['EVEX_encodings']
for item in objects:
    assert pin((R / item['path']).read_bytes()) == {k: item[k] for k in ('bytes', 'sha256')}
    inputs[item['path']] = {k: item[k] for k in ('bytes', 'sha256')}
new_text = {}
name = 'src/gallium/targets/lavapipe/meson.build'
text = (S / name).read_text()
for item in old['native_objects']:
    assert text.count(str(R / item['path'])) == 1
# Keep D3/pureF and replace the three S row objects and old adapter with seven
# actual new CPP/layer/adapter objects. Every path remains an external object.
lines = text.splitlines(True)
indices = [i for i, line in enumerate(lines) if any(str(R / item['path']) in line for item in old['native_objects'][4:])]
assert len(indices) == 4 and indices == list(range(indices[0], indices[0]+4))
indent = lines[indices[0]][:len(lines[indices[0]])-len(lines[indices[0]].lstrip())]
replacement = [indent + "'" + str(R / item['path']) + "',\n" for item in objects[4:]]
lines[indices[0]:indices[-1]+1] = replacement
new_text[name] = ''.join(lines)
name = 'src/gallium/drivers/llvmpipe/lp_state_cs.c'
text = (S / name).read_text()
changes = [
    ('   bool native_simd_packets;\n', ''),
    ('''/* Pure policy only: retain serial-pool packets; expose one row per worker
 * when all Transport rows fit the existing pool. No pool setting changes.
 */
static bool ATTRIBUTE_NOINLINE
sirius_native_transport_simd_policy(unsigned rows, unsigned workers)
{
   return workers == 0 || rows > workers;
}
''', ''),
    ('''   if (job_info->native_simd_packets) {
      /* Full packets cover eight rows; each remainder row owns one task. */
      unsigned full_packets = job_info->grid_size[0] / 8;
      grid_z = grid_y = 0;
      grid_x = (unsigned)iter_idx < full_packets ? (unsigned)iter_idx * 8 :
         full_packets * 8 + ((unsigned)iter_idx - full_packets);
   } else if (job_info->use_iters) {''', '   if (job_info->use_iters) {'),
    ('''         job_info.native_simd_packets = !native_zero_active && sirius_native_simd_supported() &&
            (job_info.native_cpu_stage != 1 ||
             sirius_native_transport_simd_policy(job_info.grid_size[0], screen->cs_tpool->num_threads));
''', ''),
    ('''   /* Native full packets and individual tails retain original workgroup IDs.
    * Other CPUs use the original scalar row tasks; ordinary JIT is unchanged. */''',
     '''   /* SIMD now evaluates independent nodes within each native row. Keep
    * one original logical row per task on every supported CPU and pool. */'''),
    ('''   } else if (job_info.native_simd_packets) {
      scheduled_tasks = job_info.grid_size[0] / 8 + job_info.grid_size[0] % 8;
''', '')]
for a, b in changes:
    assert text.count(a) == 1, a
    text = text.replace(a, b, 1)
assert 'native_simd_packets' not in text and 'sirius_native_transport_simd_policy' not in text
new_text[name] = text
name = 'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h'
text = (S / name).read_text()
a = '''   if (job->native_simd_packets)
      groups = x < (job->grid_size[0] / 8) * 8 ? 8 : 1;
'''
assert text.count(a) == 1
new_text[name] = text.replace(a, '', 1)
for name, a, b in (
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8S', 'SIR-CPU-25.2.8T'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-zero-term-v17', 'Sirius-native-cpu-field-simd-v18')):
    text = (S / name).read_text(); assert text.count(a) == 1
    new_text[name] = text.replace(a, b, 1)
patch, after, modified = '', dict(before), {}
for name, text in new_text.items():
    original = (S / name).read_text()
    patch += ''.join(difflib.unified_diff(original.splitlines(True), text.splitlines(True), n=3,
        fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(text.encode()); modified[name] = {'basis': before[name], 'candidate': after[name]}
for path in (W / 'prepare_provider_patch.py', W / 'executor-derivation.json', W / 'schedule-preparation.json',
             W / 'object-code-inspection.json', R / 'out/native-cpu-field-compatibility/result.json',
             R / 'out/native-cpu-field-compatibility/run/owner.json',
             R / 'out/native-cpu-zero-term/provider/source-contract.json'):
    inputs[str(path.relative_to(R))] = pin(path.read_bytes())
assert doc(R / 'out/native-cpu-field-compatibility/result.json')['passed']
contract = {'basis_source_identity': old['full_candidate_source_identity'],
    'full_candidate_source_identity': canonical(after), 'modified_files': modified,
    'native_objects': objects, 'source_inputs': inputs, 'patch': pin(patch.encode()),
    'new_CPP_compilations': 4, 'new_ISPC_compilations': 3, 'reused_native_objects': 4,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'scope': 'Within-row field SIMD with original scalar control and D fallback, actual OS/ISA, pureF/B admission, all9380instructions. Original logical row tasks/count64 retained, old row-SIMD scheduling retired; no threshold tuning, FP/cache/governor/host/SPV changes. Full38 Transport exact comparison passed. Build/source only; C8/all32, ABBA22+Decimal/all8 and original frame remain required.'}
P.mkdir(exist_ok=True)
(P / 'provider.patch').write_text(patch)
(P / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(path.relative_to(S)): pin(path.read_bytes()) for path in S.rglob('*') if path.is_file()}
print(json.dumps({'modified_provider_files': len(modified), 'native_objects': len(objects), 'shared_Mesa_modified': False}))
