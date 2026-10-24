"""Link the exact same shaders under the admitted 32-bit offset model."""
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
original = json.loads((root / 'out/native-cpu-simd-rows/provider/source-contract.json').read_text())
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

new_text = {}
name = 'src/gallium/targets/lavapipe/meson.build'
text = (source_root / name).read_text()
for old, new in zip(original['native_objects'][4:], objects[4:]):
    a, b = str(root / old['path']), str(root / new['path'])
    assert text.count(a) == 1
    text = text.replace(a, b)
new_text[name] = text
for name, a, b in [
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8G', 'SIR-CPU-25.2.8H'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-simd-rows-v7', 'Sirius-native-cpu-simd-address32-v8'),
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
inputs = {**simd['source_inputs'], **adapter['source_inputs']}
contract = {'scope': 'Exact G uniform shader source/all9380calls/all5programs under documented32-bit offsets with64-bit pointers. Full B capacity bounds keep every admitted view below2^31. All other ISA/FP/math flags, group/query counts, pure admission, SIMD8/scalar tails, fallback and FP scope unchanged. No source arithmetic rewrite, eligibility reduction or numerical/benefit/frame/science/release claim.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()), 'source_inputs': inputs,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp'}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'changed_provider_files': len(modified), 'original_D_F_objects_reused': 4}))
