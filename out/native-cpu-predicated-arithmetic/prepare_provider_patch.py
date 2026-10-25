"""Replace only the three raw objects and private cache/provider identities."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
V = R / 'out/native-cpu-product-inline/provider'
doc = lambda p: json.loads(p.read_text())
pin = lambda data: {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
canonical = lambda value: pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())
assert not (P / 'provider.patch').exists()
old = doc(V / 'source-contract.json')
before = {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
assert canonical(before) == old['full_candidate_source_identity']
compiled = doc(W / 'ispc-predicated-compile-result.json')
assert compiled['all_objects_compiled'] and compiled['new_compilations'] == 3
assert doc(W / 'predicated-code-inspection.json')['passed']
objects = old['native_objects'][:4] + [row['object'] for row in compiled['objects']] + old['native_objects'][7:]
assert len(objects) == 8
for row in objects:
    assert pin((R / row['path']).read_bytes()) == {k: row[k] for k in ('bytes', 'sha256')}
texts = {}
name = 'src/gallium/targets/lavapipe/meson.build'; text = (S / name).read_text()
for a, b in zip(old['native_objects'][4:7], objects[4:7]):
    assert text.count(str(R / a['path'])) == 1
    text = text.replace(str(R / a['path']), str(R / b['path']), 1)
texts[name] = text
for name, a, b in (
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-product-inline-v21', 'Sirius-native-cpu-predicated-arithmetic-v22'),
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8W', 'SIR-CPU-25.2.8X')):
    text = (S / name).read_text(); assert text.count(a) == 1; texts[name] = text.replace(a, b, 1)
patch = ''; after = dict(before); modified = {}
for name, text in texts.items():
    original = (S / name).read_text()
    patch += ''.join(difflib.unified_diff(original.splitlines(True), text.splitlines(True), n=3,
        fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(text.encode()); modified[name] = {'basis': before[name], 'candidate': after[name]}
inputs = {**old['source_inputs'], **compiled['source_inputs']}
for path in (V / 'source-contract.json', V / 'provider-build-result.json', V / 'linked-native-proof.json',
             W / 'ispc-predicated-compile-result.json', W / 'compile/ispc-predicated/owner.json',
             W / 'object-code-inspection-predicated.json', W / 'predicated-code-inspection.json', W / 'prepare_provider_patch.py'):
    inputs[str(path.relative_to(R))] = pin(path.read_bytes())
for name, expected in inputs.items():
    assert pin((R / name).read_bytes()) == expected, name
contract = {'basis_source_identity': old['full_candidate_source_identity'],
    'full_candidate_source_identity': canonical(after), 'modified_files': modified, 'patch': pin(patch.encode()),
    'native_objects': objects, 'source_inputs': inputs, 'new_native_compilations': 3,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'scope': 'X replaces only W raw3 with three pure arithmetic primitives using explicit lane result selection and the original ordered expressions. U scalar3, pureF, Iadapter, Q scheduler, V deferred NIR, full canonical program/SPIRV and all admission/fallback/FP/count/field/rollback/rounding contracts remain. Only meson object paths and private cache/provider IDs change. Additional unselected-lane FP status remains inside the existing masked/restored scope. Actual source/build/code only; fresh full38, original C8, sole ABBA22+Decimal/all8/frame remain required; no science/release/main claim.'}
P.mkdir(exist_ok=True)
(P / 'provider.patch').write_text(patch)
(P / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_provider_files': 3, 'new_raw_objects': 3, 'shared_Mesa_modified': False}))
