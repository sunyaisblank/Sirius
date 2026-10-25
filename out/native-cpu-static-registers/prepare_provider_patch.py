"""Replace only the three raw objects and private cache/provider identities."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
V = R / 'out/native-cpu-direct-scratch/provider'
doc = lambda p: json.loads(p.read_text())
pin = lambda data: {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
canonical = lambda value: pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())
assert not (P / 'provider.patch').exists()
old = doc(V / 'source-contract.json')
before = {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
assert canonical(before) == old['full_candidate_source_identity']
compiled = doc(W / 'ispc-registers-compile-result.json')
assert compiled['all_objects_compiled'] and compiled['new_compilations'] == 3
assert doc(W / 'static-register-code-inspection.json')['passed']
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
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-direct-scratch-v24', 'Sirius-native-cpu-static-registers-v25'),
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8Z', 'SIR-CPU-25.2.8AA')):
    text = (S / name).read_text(); assert text.count(a) == 1; texts[name] = text.replace(a, b, 1)
patch = ''; after = dict(before); modified = {}
for name, text in texts.items():
    original = (S / name).read_text()
    patch += ''.join(difflib.unified_diff(original.splitlines(True), text.splitlines(True), n=3,
        fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(text.encode()); modified[name] = {'basis': before[name], 'candidate': after[name]}
inputs = {**old['source_inputs'], **compiled['source_inputs']}
for path in (V / 'source-contract.json', V / 'provider-build-result.json', V / 'linked-native-proof.json',
             W / 'ispc-registers-compile-result.json', W / 'compile/ispc-registers/owner.json',
             W / 'object-code-inspection-registers.json', W / 'static-register-code-inspection.json', W / 'prepare_provider_patch.py'):
    inputs[str(path.relative_to(R))] = pin(path.read_bytes())
for name, expected in inputs.items():
    assert pin((R / name).read_bytes()) == expected, name
contract = {'basis_source_identity': old['full_candidate_source_identity'],
    'full_candidate_source_identity': canonical(after), 'modified_files': modified, 'patch': pin(patch.encode()),
    'native_objects': objects, 'source_inputs': inputs, 'new_native_compilations': 3,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'scope': 'AA replaces only Z raw3 with compile-time uniform canonical register metadata. Exact whole-source inverse restores Z; all9380ordered tuples and arithmetic bodies unchanged. Actual family/register/layer guards still dominate static calls; generic fallback and per-row numerical refusal/publication retained. Same Z direct-scratch ABI, Y scheduling/FP scope, U scalar3/pureF/Iadapter, Q scheduling, V deferred NIR and original SPIRV. Only Meson object paths and private provider/cache IDs change. Source/build/code only; fresh full38/C8 and original unchanged frame required. Z broader controls are explicitly predecessor evidence; no current all8/all22 or science/release/main claim.'}
P.mkdir(exist_ok=True)
(P / 'provider.patch').write_text(patch)
(P / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_provider_files': 3, 'new_raw_objects': 3, 'shared_Mesa_modified': False}))
