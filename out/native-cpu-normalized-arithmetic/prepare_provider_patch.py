"""Reuse the raw objects; replace the combined arithmetic adapter and private identities."""
from pathlib import Path
import difflib
import hashlib
import json

R = Path.cwd(); W = Path(__file__).resolve().parent; P = W / 'provider'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
V = R / 'out/native-cpu-arithmetic-closure/provider'
doc = lambda p: json.loads(p.read_text())
pin = lambda data: {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
canonical = lambda value: pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())
assert not (P / 'provider.patch').exists()
old = doc(V / 'source-contract.json')
before = {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
assert canonical(before) == old['full_candidate_source_identity']
compiled=doc(R/'out/native-cpu-arithmetic-closure/ispc-dense-compile-result.json')
assert compiled['all_objects_compiled'] and compiled['new_compilations']==3
cpp=doc(W/'cpp-build-result.json');assert cpp['passed'] and cpp['new_relocatable_links']==1
assert doc(W/'cpp-code-inspection.json')['passed']
assert doc(W/'domain/run/owner.json')['passed'] and doc(W/'domain/run/stdout.log')['compared_words']==17781120
assert cpp['adapter_input_code_association']['passed']
objects=old['native_objects'][:7]+[cpp['merged_adapter']]
assert objects[4:7]==[v['object'] for v in compiled['objects']]
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
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-arithmetic-closure-v31', 'Sirius-native-cpu-normalized-arithmetic-v32'),
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8AG', 'SIR-CPU-25.2.8AH')):
    text = (S / name).read_text(); assert text.count(a) == 1; texts[name] = text.replace(a, b, 1)
patch = ''; after = dict(before); modified = {}
for name, text in texts.items():
    original = (S / name).read_text()
    patch += ''.join(difflib.unified_diff(original.splitlines(True), text.splitlines(True), n=3,
        fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(text.encode()); modified[name] = {'basis': before[name], 'candidate': after[name]}
inputs = {**old['source_inputs'], **doc(W/'prepared-cpp.json')['inputs'],**doc(W/'prepared-domain.json')['inputs']}
for path in (V/'source-contract.json',V/'provider-build-result.json',V/'linked-native-proof.json',W/'cpp-build-result.json',W/'cpp/build/owner.json',W/'cpp-code-inspection.json',W/'domain-build-result.json',W/'domain/run/owner.json',W/'domain/run/stdout.log',W/'prepare_provider_patch.py'):
    inputs[str(path.relative_to(R))] = pin(path.read_bytes())
for name, expected in inputs.items():
    assert pin((R / name).read_bytes()) == expected, name
contract = {'basis_source_identity': old['full_candidate_source_identity'],
    'full_candidate_source_identity': canonical(after), 'modified_files': modified, 'patch': pin(patch.encode()),
    'native_objects': objects, 'source_inputs': inputs, 'new_native_compilations':1,'new_ISPC_compilations':0,'reused_raw_objects':3,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'scope': 'AH reuses exact AG raw3 and all9380tuples. One new strict-FP C++ helper adds guarded normalized fast domains, retaining the complete original fallback for every other input. Public final2^120 error/normalization/refusal and original identity/order remain; reduced private predicates follow from explicit spacing bounds, with Multiply highs<=2^54 implying rounded low products<=2^62. Full-mask stores avoid old-output reads and normal exponents avoid unused reconstruction. Original I/new CPP relocatable code join passed, plus two scoped original-scalar-reference fixtures; extended fixture17781120words includes full8-lane normalized/mixed/fallback/masks/aliases. FP/precision/eligibility/scheduler/queue/layout/budgets remain. Current host separately reuses20already-owned lower increments without new shader requests. Fresh linked targets/full38/C8/affected original control5/frame required; no performance/science/release claim.'}
P.mkdir(exist_ok=True)
(P / 'provider.patch').write_text(patch)
(P / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(p.relative_to(S)): pin(p.read_bytes()) for p in S.rglob('*') if p.is_file()}
print(json.dumps({'modified_provider_files': 3, 'new_raw_objects': 0, 'reused_raw_objects': 3, 'shared_Mesa_modified': False}))
