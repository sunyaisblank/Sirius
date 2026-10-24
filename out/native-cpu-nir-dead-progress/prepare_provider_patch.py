"""Compose one NIR DCE pass and a new private identity; reuse all L objects."""
from pathlib import Path
import difflib
import hashlib
import json
import runpy

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
previous = root / 'out/native-cpu-simd-endpoint-regions'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(value):
    return pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


assert not (provider / 'provider.patch').exists()
original = json.loads((previous / 'provider/source-contract.json').read_text())
before = {str(path.relative_to(source_root)): pin(path.read_bytes())
          for path in source_root.rglob('*') if path.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
objects = original['native_objects']
assert len(objects) == 8
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {key: item[key] for key in ('bytes', 'sha256')}
simd = json.loads((previous / 'ispc-uniform-compile-result.json').read_text())
owner = json.loads((previous / 'compile/ispc-uniform/owner.json').read_text())
assert simd['all_objects_compiled'] and owner['passed'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes']
assert [item['object'] for item in simd['objects']] == objects[4:7]
proof = json.loads((work / 'source-proof.json').read_text())
assert proof['native_objects_reused_without_recompilation'] == objects
assert proof['new_native_objects'] == 0 and not proof['patch_applied']
namespace = runpy.run_path(str(work / 'prepare_patch.py'))
nir_name = proof['changed_file']
assert pin((source_root / nir_name).read_bytes()) == proof['basis']
assert pin(namespace['new'].encode()) == proof['candidate']
new_text = {nir_name: namespace['new']}
for name, old, new in [
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8L', 'SIR-CPU-25.2.8M'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-simd-endpoint-regions-v12',
     'Sirius-native-cpu-nir-dead-progress-v13'),
]:
    text = (source_root / name).read_text()
    assert text.count(old) == 1
    new_text[name] = text.replace(old, new, 1)
assert len(new_text) == 3

patch, after, modified = '', dict(before), {}
for name, new in new_text.items():
    old = (source_root / name).read_text()
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=3,
                                         fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(new.encode())
    modified[name] = {'basis': before[name], 'candidate': after[name]}
inputs = dict(original['source_inputs'])
source_prefix = str(source_root.relative_to(root)) + '/'
for name, value in proof['source_inputs'].items():
    # The candidate full-tree identity seals this changed file after apply.
    # Do not retain its old basis hash as a future runtime-input expectation.
    if name == source_prefix + nir_name:
        continue
    assert pin((root / name).read_bytes()) == value, name
    inputs[name] = value
inputs.pop(source_prefix + nir_name, None)
for name in ('prepare_patch.py', 'provider.patch', 'source-proof.json'):
    p = work / name
    inputs[str(p.relative_to(root))] = pin(p.read_bytes())
for name in ('provider/source-contract.json', 'ispc-uniform-compile-result.json',
             'compile/ispc-uniform/owner.json'):
    p = previous / name
    inputs[str(p.relative_to(root))] = pin(p.read_bytes())
assert all(pin((root / name).read_bytes()) == value for name, value in inputs.items())
contract = {
    'scope': 'M adds one end-of-main-loop nir_opt_dce cleanup in lp_build_opt_nir to eliminate dead SSA rewrite cones. All8 L native objects reused at their exact paths/bytes; no ISPC/C++ regeneration or compilation. L Endpoint regions, K Transport/Dopri, I guarded SIMD/scalar-tail policy, separate zero-header no-op, full positive B admission, original FP/layout/counts/limits/eligibility and ordinary fallback remain unchanged. Actual prior samples identify NIR passes; precise module/cycle and timeout cure remain unproven. No numerical/benefit/frame/scientific/release claim.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()), 'source_inputs': inputs,
    'new_native_compilations': 0, 'adapter_recompiled': False,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'nir_dead_progress_source_proof': pin((work / 'source-proof.json').read_bytes()),
}
provider.mkdir(exist_ok=True)
(provider / 'provider.patch').write_text(patch)
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
reuse = {'scope': 'Truthful reuse receipt; not a compiler result.',
         'basis_provider_contract': pin((previous / 'provider/source-contract.json').read_bytes()),
         'native_objects': objects, 'all8_object_bytes_verified': True,
         'new_native_compilations': 0, 'new_native_object_links_or_copies': 0}
assert not (work / 'reused-native-objects.json').exists()
(work / 'reused-native-objects.json').write_text(json.dumps(reuse, indent=2) + '\n')
assert before == {str(path.relative_to(source_root)): pin(path.read_bytes())
                  for path in source_root.rglob('*') if path.is_file()}
print(json.dumps({'patch_written': True, 'changed_provider_files': 3,
                  'native_objects_reused': 8, 'new_native_compilations': 0,
                  'shared_source_modified': False}))
