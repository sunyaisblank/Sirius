"""Link all three actual S objects into exact R Mesa; retain Q scheduling."""
from pathlib import Path
import difflib
import hashlib
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
prior = root / 'out/native-cpu-accumulate-inline'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'


def doc(path):
    return json.loads(path.read_text())


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(value):
    return pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


assert not (provider / 'provider.patch').exists()
original = doc(prior / 'provider/source-contract.json')
before = {str(path.relative_to(source_root)): pin(path.read_bytes())
          for path in source_root.rglob('*') if path.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
compiled = doc(work / 'ispc-zero-compile-result.json')
owner = doc(work / 'compile/ispc-zero/owner.json')
inline = doc(work / 'zero-code-inspection.json')
isa = doc(work / 'object-code-inspection-zero.json')
assert compiled['all_objects_compiled'] and compiled['all_prepared_source_bytes_unchanged']
assert compiled['new_compilations'] == 3 and compiled['reused_compiled_objects'] == 0
assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
assert inline['passed'] and inline['compile_result'] == pin((work / 'ispc-zero-compile-result.json').read_bytes())
assert [row['stage'] for row in compiled['objects']] == ['transport', 'endpoint', 'dopri_phase']
objects = original['native_objects'][:4] + [row['object'] for row in compiled['objects']] + original['native_objects'][7:]
assert len(objects) == 8
for row in compiled['objects']:
    slug = row['stage']
    assert row['argv_executed'] and row['exit'] == 0
    assert inline['stages'][slug]['all_RTAccumulateTerm_ABI_calls_removed']
    assert inline['stages'][slug]['candidate']['object'] == isa['objects'][slug]['object'] == row['object']
    assert not isa['objects'][slug]['forbidden_fma_fp_control_x87_found']
    assert not isa['objects'][slug]['forbidden_non_v3_extensions'] and isa['objects'][slug]['EVEX_encodings'] == 0
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {key: item[key] for key in ('bytes', 'sha256')}
meson = 'src/gallium/targets/lavapipe/meson.build'
text = (source_root / meson).read_text()
for old, new in zip(original['native_objects'][4:7], objects[4:7]):
    a, b = str(root / old['path']), str(root / new['path'])
    assert text.count(a) == 1
    text = text.replace(a, b, 1)
new_text = {meson: text}
for name, old, new in (
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8R', 'SIR-CPU-25.2.8S'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-accumulate-inline-v16',
     'Sirius-native-cpu-zero-term-v17'),
):
    old_text = (source_root / name).read_text()
    assert old_text.count(old) == 1
    new_text[name] = old_text.replace(old, new, 1)
patch, after, modified = '', dict(before), {}
for name, new in new_text.items():
    old = (source_root / name).read_text()
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=3,
        fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(new.encode())
    modified[name] = {'basis': before[name], 'candidate': after[name]}
inputs = {**original['source_inputs'], **compiled['source_inputs']}
for name in new_text:
    inputs.pop(str((source_root / name).relative_to(root)), None)
for name in ('source-proof.json', 'zero-term.patch', 'ispc-derivation.json',
    'ispc-zero-compile-result.json', 'compile/ispc-zero/owner.json',
    'zero-code-inspection.json', 'object-code-inspection-zero.json', 'prepare_provider_patch.py'):
    path = work / name
    inputs[str(path.relative_to(root))] = pin(path.read_bytes())
for name in ('provider/source-contract.json', 'provider/linked-native-proof.json'):
    path = prior / name
    inputs[str(path.relative_to(root))] = pin(path.read_bytes())
assert all(pin((root / name).read_bytes()) == expected for name, expected in inputs.items())
contract = {'scope': 'S links three genuinely compiled complete ISPC objects with only the reviewed exact signed-zero arithmetic arm added; five R scalar/admission/adapter objects reused directly. Q scheduler/ISA/zero-active/DCE/guards/FP/body9380calls/counts/governors remain unchanged. P host/frame products reused. New Transport raw code requires fresh full38 compatibility before C8; no historical identical-Transport execution shortcut and no numerical/benefit/frame/scientific/release claim.',
    'basis_source_identity': original['full_candidate_source_identity'], 'modified_files': modified,
    'full_candidate_source_identity': canonical(after), 'native_objects': objects,
    'patch': pin(patch.encode()), 'source_inputs': inputs, 'new_native_compilations': 3,
    'new_CPP_compilations': 0, 'adapter_recompiled': False,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'zero_term_source_proof': pin((work / 'source-proof.json').read_bytes()),
    'zero_term_actual_code_proof': pin((work / 'zero-code-inspection.json').read_bytes()),
    'required_before_C8': 'out/native-cpu-zero-term-compatibility/result.json passed for these exact S objects'}
provider.mkdir(exist_ok=True)
(provider / 'provider.patch').write_text(patch)
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(path.relative_to(source_root)): pin(path.read_bytes())
                  for path in source_root.rglob('*') if path.is_file()}
print(json.dumps({'changed_provider_files': 3, 'new_actual_S_objects': 3,
                  'direct_R_object_reuse': 5, 'shared_Mesa_modified': False, 'executed_shader': False}))
