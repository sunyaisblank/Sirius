"""Prepare CPU-gated full SIMD packets and independent scalar tail tasks."""
from pathlib import Path
import difflib
import hashlib
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
provider.mkdir(exist_ok=True)


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(value):
    return pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after)


assert not (provider / 'provider.patch').exists()
original = json.loads((root / 'out/native-cpu-simd-address32/provider/source-contract.json').read_text())
before = {str(p.relative_to(source_root)): pin(p.read_bytes())
          for p in source_root.rglob('*') if p.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
simd = json.loads((work / 'ispc-uniform-compile-result.json').read_text())
adapter = json.loads((work / 'adapter-compile-result.json').read_text())
assert simd['all_objects_compiled'] and adapter['all_prepared_source_bytes_unchanged']
for name in ('ispc-uniform', 'adapter'):
    owner = json.loads((work / 'compile' / name / 'owner.json').read_text())
    assert owner['passed'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes']
objects = original['native_objects'][:4] + [v['object'] for v in simd['objects']] + [adapter['objects'][0]['object']]
assert len(objects) == 8
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {k: item[k] for k in ('bytes', 'sha256')}

cs_name = 'src/gallium/drivers/llvmpipe/lp_state_cs.c'
cs = replace((source_root / cs_name).read_text(),
    '   bool native_preflight_passed;\n',
    '   bool native_preflight_passed;\n   bool native_simd_packets;\n')
cs = replace(cs, '   if (job_info->use_iters) {', '''   if (job_info->native_simd_packets) {
      /* Full packets cover eight rows; each remainder row owns one task. */
      unsigned full_packets = job_info->grid_size[0] / 8;
      grid_z = grid_y = 0;
      grid_x = (unsigned)iter_idx < full_packets ? (unsigned)iter_idx * 8 :
         full_packets * 8 + ((unsigned)iter_idx - full_packets);
   } else if (job_info->use_iters) {''')
cs = replace(cs, '         job_info.native_preflight_passed = true;',
    '         job_info.native_preflight_passed = true;\n'
    '         job_info.native_simd_packets = sirius_native_simd_supported();')
cs = replace(cs, '''   /* Native SIMD packets reuse the original packed task IDs after pure
    * admission. Original grid, query counts and fallback key stay intact. */
   if (intergroup || job_info.native_preflight_passed) {''', '''   /* Native full packets and individual tails retain original workgroup IDs.
    * Other CPUs use the original scalar row tasks; ordinary JIT is unchanged. */
   if (job_info.native_simd_packets) {
      scheduled_tasks = job_info.grid_size[0] / 8 + job_info.grid_size[0] % 8;
   } else if (intergroup) {''')
bridge_name = 'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h'
bridge = replace((source_root / bridge_name).read_text(),
    '#include "sirius_native_bridge.h"', '''#include "sirius_native_bridge.h"
#include "util/u_cpu_detect.h"

/* Mesa checks AVX, OSXSAVE and enabled XMM/YMM state. The compiler's required
 * core-avx2 device also uses v3 integer extensions; guard those here too.
 * Every other CPU keeps the original scalar path, including direct callers.
 */
bool ATTRIBUTE_NOINLINE
sirius_native_simd_supported(void)
{
   return util_get_cpu_caps()->has_avx2 && __builtin_cpu_supports("x86-64-v3");
}''')
bridge = replace(bridge,
    '   uint32_t groups = job->intergroup_row_owner ? MIN2(8, job->grid_size[0] - x) : 1;',
    '''   uint32_t groups = job->intergroup_row_owner ? MIN2(8, job->grid_size[0] - x) : 1;
   if (job->native_simd_packets)
      groups = x < (job->grid_size[0] / 8) * 8 ? 8 : 1;''')
abi_name = 'src/gallium/drivers/llvmpipe/sirius_native_bridge.h'
abi = replace((source_root / abi_name).read_text(), '\n#undef SIRIUS_NATIVE_PACKET_ARGUMENTS',
    '\n// Pure runtime ISA/OS-state query; scalar baseline remains supported.\n'
    'bool sirius_native_simd_supported(void) SIRIUS_NATIVE_NOEXCEPT;\n'
    '\n#undef SIRIUS_NATIVE_PACKET_ARGUMENTS')
meson_name = 'src/gallium/targets/lavapipe/meson.build'
meson = (source_root / meson_name).read_text()
for old, new in zip(original['native_objects'][4:], objects[4:]):
    meson = replace(meson, str(root / old['path']), str(root / new['path']))
new_text = {cs_name: cs, bridge_name: bridge, abi_name: abi, meson_name: meson}
for name, a, b in [
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8H', 'SIR-CPU-25.2.8I'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-simd-address32-v8', 'Sirius-native-cpu-simd-packet-policy-v9'),
]:
    new_text[name] = replace((source_root / name).read_text(), a, b)
patch = ''
after = dict(before)
modified = {}
for name, new in new_text.items():
    old = (source_root / name).read_text()
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=0,
                                         fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(new.encode())
    modified[name] = {'basis': before[name], 'candidate': after[name]}
(provider / 'provider.patch').write_text(patch)
inputs = {**simd['source_inputs'], **adapter['source_inputs']}
contract = {'scope': 'Optional core-avx2 full8 packets guarded by Mesa OS-enabled AVX2 and x86-64-v3 CPU level, exactG/H source/all9380calls/all5programs, original F scalar1 on every otherx86CPU. Each tail row remains a separate parallel scalar task. Pure admission, FP guard, all original grid/query/counts/fallback, 32-bit bounded offsets and64-bit pointers unchanged. Direct C ABI guarded too; no FMA/fast math despite hardware support. No numerical/benefit/frame/scientific/release claim.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()), 'source_inputs': inputs,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp'}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'changed_provider_files': len(modified), 'original_D_F_objects_reused': 4}))
