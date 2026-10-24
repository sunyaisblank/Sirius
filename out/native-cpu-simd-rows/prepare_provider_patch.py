"""Prepare eight-row native scheduling and typed entrypoints; do not apply/build."""
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


def canonical(records):
    return pin(json.dumps(records, sort_keys=True, separators=(',', ':')).encode())


def replace(text, before, after):
    assert text.count(before) == 1, before
    return text.replace(before, after)


assert not (provider / 'provider.patch').exists()
original = json.loads((root / 'out/native-cpu-dispatch-admission/provider/source-contract.json').read_text())
before = {str(p.relative_to(source_root)): pin(p.read_bytes())
          for p in source_root.rglob('*') if p.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
simd = json.loads((work / 'ispc-uniform-compile-result.json').read_text())
adapter = json.loads((work / 'adapter-compile-result.json').read_text())
assert simd['all_objects_compiled'] and adapter['all_prepared_source_bytes_unchanged']
for name in ('ispc-uniform', 'adapter'):
    owner = json.loads((work / 'compile' / name / 'owner.json').read_text())
    assert owner['passed'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes']
assert adapter['objects'][0]['exit'] == 0
objects = original['native_objects'] + [item['object'] for item in simd['objects']] + [adapter['objects'][0]['object']]
assert len(objects) == 8
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {k: item[k] for k in ('bytes', 'sha256')}

cs_name = 'src/gallium/drivers/llvmpipe/lp_state_cs.c'
cs = replace((source_root / cs_name).read_text(),
    '   if (intergroup) {\n      job_info.intergroup_row_owner = true;',
    '   /* Native SIMD packets reuse the original packed task IDs after pure\n'
    '    * admission. Original grid, query counts and fallback key stay intact. */\n'
    '   if (intergroup || job_info.native_preflight_passed) {\n      job_info.intergroup_row_owner = true;')
bridge_name = 'src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h'
bridge = (source_root / bridge_name).read_text()
for slug in ('transport', 'endpoint', 'dopri_phase'):
    bridge = replace(bridge, f'sirius_native_{slug}_try_packet(', f'sirius_native_{slug}_rows_try_packet(')
abi_name = 'src/gallium/drivers/llvmpipe/sirius_native_bridge.h'
exports = '\n// Full8 SIMD packets; original scalar execution for every smaller tail.\n'
for slug in ('transport', 'endpoint', 'dopri_phase'):
    exports += f'bool sirius_native_{slug}_rows_try_packet(SIRIUS_NATIVE_PACKET_ARGUMENTS) SIRIUS_NATIVE_NOEXCEPT;\n'
abi = replace((source_root / abi_name).read_text(), '\n#undef SIRIUS_NATIVE_PACKET_ARGUMENTS', exports + '\n#undef SIRIUS_NATIVE_PACKET_ARGUMENTS')
meson_name = 'src/gallium/targets/lavapipe/meson.build'
old_last = "    '" + str(root / objects[3]['path']) + "'\n"
new_last = old_last.rstrip('\n') + ',\n' + ',\n'.join("    '" + str(root / item['path']) + "'" for item in objects[4:]) + '\n'
meson = replace((source_root / meson_name).read_text(), old_last, new_last)
new_text = {cs_name: cs, bridge_name: bridge, abi_name: abi, meson_name: meson}
for name, a, b in [
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8F', 'SIR-CPU-25.2.8G'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-dispatch-admission-v6', 'Sirius-native-cpu-simd-rows-v7'),
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
for name in ('ispc-uniform-compile-result.json', 'adapter-compile-result.json', 'object-code-inspection-uniform.json'):
    p = work / name
    inputs[str(p.relative_to(root))] = pin(p.read_bytes())
contract = {'scope': 'Private eight-row SIMD native packets with all5canonicalprograms/all9380originalcalls. Same full F/B admission, original scalar1..7 tails, original group/query counts, fallback intergroup key/JIT, pending variant, FP guard and complete T/E/DP contracts. No public source adoption or scientific/frame/benefit claim.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()), 'source_inputs': inputs,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'limits': 'Source/build preparation only. Actual all32 readbacks, original soleABBA/all22/all8/frame and full science/release remain pending.'}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'changed_provider_files': len(modified), 'reused_original_objects': 4, 'new_simd_and_adapter_objects': 4}))
