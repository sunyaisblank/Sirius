"""Link uniform structural indices using the unchanged I packet policy/adapter."""
from pathlib import Path
import difflib
import hashlib
import json
import runpy

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
provider.mkdir(exist_ok=True)


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(value):
    return pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


assert not (provider / 'provider.patch').exists()
previous = root / 'out/native-cpu-simd-structural-uniform'
adapter_basis = root / 'out/native-cpu-simd-packet-policy'
original = json.loads((previous / 'provider/source-contract.json').read_text())
before = {str(p.relative_to(source_root)): pin(p.read_bytes())
          for p in source_root.rglob('*') if p.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
simd = json.loads((work / 'ispc-uniform-compile-result.json').read_text())
owner = json.loads((work / 'compile/ispc-uniform/owner.json').read_text())
assert simd['all_objects_compiled'] and owner['passed'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes']
adapter = json.loads((adapter_basis / 'adapter-compile-result.json').read_text())
owner = json.loads((adapter_basis / 'compile/adapter/owner.json').read_text())
assert adapter['all_prepared_source_bytes_unchanged'] and owner['passed'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes']
assert adapter['objects'][0]['object'] == original['native_objects'][7]
# Only opaque word pointers/counts cross this ABI; private context layout is
# entirely inside each new ISPC object. The guarded I adapter is unchanged.
objects = original['native_objects'][:4] + [v['object'] for v in simd['objects']] + [original['native_objects'][7]]
assert len(objects) == 8
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {k: item[k] for k in ('bytes', 'sha256')}

name = 'src/gallium/targets/lavapipe/meson.build'
text = (source_root / name).read_text()
for old, new in zip(original['native_objects'][4:7], objects[4:7]):
    a, b = str(root / old['path']), str(root / new['path'])
    assert text.count(a) == 1
    text = text.replace(a, b)
new_text = {name: text}
for name, a, b in [
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8K', 'SIR-CPU-25.2.8L'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-simd-structural-uniform-v11', 'Sirius-native-cpu-simd-endpoint-regions-v12'),
]:
    text = (source_root / name).read_text()
    assert text.count(a) == 1
    new_text[name] = text.replace(a, b)
# Compose the separately proved zero-header no-op into this same private
# candidate, so its source is qualified together with Endpoint regions.
zero = root / 'out/native-cpu-zero-active-prepare'
zero_namespace = runpy.run_path(str(zero / 'prepare_patch.py'))
zero_proof = json.loads((zero / 'source-proof.json').read_text())
assert zero_proof['positive_B_guard_unchanged']
assert zero_proof['zero_header_classifier_separate_from_positive_admission']
for name, text in zero_namespace['new_sources'].items():
    relative = 'src/gallium/drivers/llvmpipe/' + name
    assert relative not in new_text
    assert pin((source_root / relative).read_bytes()) == zero_proof['changed_files'][relative]['basis']
    assert pin(text.encode()) == zero_proof['changed_files'][relative]['candidate']
    new_text[relative] = text
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
inputs = {**adapter['source_inputs'], **simd['source_inputs']}
for name in ('adapter-compile-result.json', 'compile/adapter/owner.json', 'rows_adapter.cpp', 'rows_bridge.h'):
    p = adapter_basis / name
    inputs[str(p.relative_to(root))] = pin(p.read_bytes())
for name in ('prepare_patch.py', 'provider.patch', 'source-proof.json'):
    path = zero / name
    inputs[str(path.relative_to(root))] = pin(path.read_bytes())
for name, value in zero_proof['source_inputs'].items():
    if not name.startswith(str(source_root.relative_to(root)) + '/') or name.removeprefix(str(source_root.relative_to(root)) + '/') not in zero_proof['changed_files']:
        assert pin((root / name).read_bytes()) == value
        inputs[name] = value
contract = {'zero_active_source_patch': {'proof': pin((zero / 'source-proof.json').read_bytes()),
    'positive_B_admission_unchanged': True, 'separate_no_op_class': True,
    'original_query_invocations_and_fallback_state_preserved': True}, 'scope': 'Fixed Endpoint induction regions with full inverse to K and exact K Transport/Dopri reuse; separate exact-zero-header no-op with original query counts/fallback state. All positive B numerical admission and arithmetic remain unchanged. All9380literalcalls/all5programs/arithmetic and real input-row identities retained. Exact I guarded adapter, Mesa CPU/v3 routing, independent scalar tails, original F baseline, full B guard/counts/governors unchanged. No numerical/benefit/frame/scientific/release claim.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()), 'source_inputs': inputs,
    'adapter_recompiled': False, 'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp'}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'changed_provider_files': len(modified), 'objects_reused': 5 + simd['reused_compiled_objects'], 'new_SIMD_objects': simd['new_compilations'], 'linked_SIMD_objects': 3}))
