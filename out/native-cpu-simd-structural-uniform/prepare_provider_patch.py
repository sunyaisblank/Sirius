"""Link uniform structural indices using the unchanged I packet policy/adapter."""
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


assert not (provider / 'provider.patch').exists()
previous = root / 'out/native-cpu-simd-output-columns'
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
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8J', 'SIR-CPU-25.2.8K'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-simd-output-columns-v10', 'Sirius-native-cpu-simd-structural-uniform-v11'),
]:
    text = (source_root / name).read_text()
    assert text.count(a) == 1
    new_text[name] = text.replace(a, b)
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
contract = {'scope': 'Uniform structural indices and contiguous private output/coefficient access, derived from accepted J with closed inverse proof. All9380literalcalls/all5programs/arithmetic and real input-row identities retained. Exact I guarded adapter, Mesa CPU/v3 routing, independent scalar tails, original F baseline, full B guard/counts/governors unchanged. No numerical/benefit/frame/scientific/release claim.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()), 'source_inputs': inputs,
    'adapter_recompiled': False, 'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp'}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'changed_provider_files': len(modified), 'objects_reused': 5, 'new_SIMD_objects': 3}))
