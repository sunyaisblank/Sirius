"""Replace only the three raw objects and private cache/provider identities."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
V = R / 'out/native-cpu-demand-submit/provider'
doc = lambda p: json.loads(p.read_text())
pin = lambda data: {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
canonical = lambda value: pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())
assert not (P / 'provider.patch').exists()
old = doc(V / 'source-contract.json')
before = {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
assert canonical(before) == old['full_candidate_source_identity']
compiled = doc(W / 'ispc-dense-compile-result.json')
assert compiled['all_objects_compiled'] and compiled['new_compilations'] == 3
assert doc(W / 'dense-code-inspection.json')['passed']
objects = old['native_objects'][:4] + [row['object'] for row in compiled['objects']] + [compiled['merged_adapter']]
assert len(objects) == 8
for row in objects:
    assert pin((R / row['path']).read_bytes()) == {k: row[k] for k in ('bytes', 'sha256')}
texts = {}
name = 'src/gallium/targets/lavapipe/meson.build'; text = (S / name).read_text()
for a, b in zip(old['native_objects'][4:], objects[4:]):
    assert text.count(str(R / a['path'])) == 1
    text = text.replace(str(R / a['path']), str(R / b['path']), 1)
texts[name] = text
for name, a, b in (
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-demand-submit-v30', 'Sirius-native-cpu-arithmetic-closure-v31'),
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8AF', 'SIR-CPU-25.2.8AG')):
    text = (S / name).read_text(); assert text.count(a) == 1; texts[name] = text.replace(a, b, 1)
patch = ''; after = dict(before); modified = {}
for name, text in texts.items():
    original = (S / name).read_text()
    patch += ''.join(difflib.unified_diff(original.splitlines(True), text.splitlines(True), n=3,
        fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(text.encode()); modified[name] = {'basis': before[name], 'candidate': after[name]}
inputs = {**old['source_inputs'], **compiled['source_inputs']}
for path in (V / 'source-contract.json', V / 'provider-build-result.json', V / 'linked-native-proof.json',
             W / 'ispc-dense-compile-result.json', W / 'compile/ispc-dense/owner.json',
             W / 'object-code-inspection-dense.json', W / 'dense-code-inspection.json', W / 'prepare_provider_patch.py', W / 'inspect_closure_objects.py'):
    inputs[str(path.relative_to(R))] = pin(path.read_bytes())
for name, expected in inputs.items():
    assert pin((R / name).read_bytes()) == expected, name
contract = {'basis_source_identity': old['full_candidate_source_identity'],
    'full_candidate_source_identity': canonical(after), 'modified_files': modified, 'patch': pin(patch.encode()),
    'native_objects': objects, 'source_inputs': inputs, 'new_native_compilations': 3,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'scope': 'AG composes the original divide/sqrt three correction steps and final residual in the same strict-FP C++ unit as checked Add/Multiply, replacing repeated bank/ABI/unpack calls with one initialized bank per compound division/sqrt. Original literal instructions and all9380tuples, row/control/refusal/rollback/publication/FP/eligibility/counts remain; original I and new CPP opcode-joined. AF on-demand queue and Y FP/scheduler retained. Same eight-object graph, three raw object paths and merged adapter plus private cache/UUID change only. Focused four-op23760call/2851200word compatibility passed; actual provider joins/full38/C8 and sole original frame required. Z broader results remain predecessor evidence, no science/release claim.'}
P.mkdir(exist_ok=True)
(P / 'provider.patch').write_text(patch)
(P / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_provider_files': 3, 'new_raw_objects': 3, 'shared_Mesa_modified': False}))
