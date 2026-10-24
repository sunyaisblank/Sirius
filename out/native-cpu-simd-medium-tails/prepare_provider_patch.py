"""Compose the medium-tail scheduler and sole genuine adapter into provider N."""
from pathlib import Path
import difflib
import hashlib
import json
import runpy

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
previous = root / 'out/native-cpu-nir-dead-progress'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(value):
    return pin(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


assert not (provider / 'provider.patch').exists()
original = json.loads((previous / 'provider/source-contract.json').read_text())
before = {str(path.relative_to(source_root)): pin(path.read_bytes())
          for path in source_root.rglob('*') if path.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
compiled = json.loads((work / 'adapter-compile-result.json').read_text())
owner = json.loads((work / 'compile/adapter/owner.json').read_text())
assert compiled['all_prepared_source_bytes_unchanged'] and owner['passed'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes']
adapter, = compiled['objects']
assert adapter['exit'] == 0
assert adapter['object']['path'] == 'bin/linux-gcc/native-cpu-simd-medium-tails/objects/rows_adapter.o'
objects = original['native_objects'][:7] + [adapter['object']]
assert len(original['native_objects']) == len(objects) == 8
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {key: item[key] for key in ('bytes', 'sha256')}
namespace = runpy.run_path(str(work / 'prepare_packet_patch.py'))
proof = json.loads((work / 'partition-proof.json').read_text())
assert proof['seven_reused_native_objects'] == objects[:7]
assert proof['adapter_only_three_group_predicates_changed'] and not proof['patch_applied']
new_text = {}
for name, text in namespace['new_sources'].items():
    relative = 'src/gallium/drivers/llvmpipe/' + name
    assert pin((source_root / relative).read_bytes()) == proof['changed_files'][relative]['basis']
    assert pin(text.encode()) == proof['changed_files'][relative]['candidate']
    new_text[relative] = text
name = 'src/gallium/targets/lavapipe/meson.build'
text = (source_root / name).read_text()
old_path = str(root / original['native_objects'][7]['path'])
new_path = str(root / adapter['object']['path'])
assert old_path != new_path and text.count(old_path) == 1
new_text[name] = text.replace(old_path, new_path, 1)
for name, old, new in [
    ('src/gallium/frontends/lavapipe/lvp_device.c', 'SIR-CPU-25.2.8M', 'SIR-CPU-25.2.8N'),
    ('src/gallium/drivers/llvmpipe/lp_screen.c', 'Sirius-native-cpu-nir-dead-progress-v13',
     'Sirius-native-cpu-simd-medium-tails-v14'),
]:
    text = (source_root / name).read_text()
    assert text.count(old) == 1
    new_text[name] = text.replace(old, new, 1)
assert len(new_text) == 5

patch, after, modified = '', dict(before), {}
for name, new in new_text.items():
    old = (source_root / name).read_text()
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=3,
                                         fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(new.encode())
    modified[name] = {'basis': before[name], 'candidate': after[name]}
inputs = {**original['source_inputs'], **compiled['source_inputs']}
source_prefix = str(source_root.relative_to(root)) + '/'
for name, value in proof['source_inputs'].items():
    if name.startswith(source_prefix) and name.removeprefix(source_prefix) in new_text:
        continue
    assert pin((root / name).read_bytes()) == value, name
    inputs[name] = value
for name in new_text:
    inputs.pop(source_prefix + name, None)
    inputs.pop(str(source_root / name), None)
for base, names in [
    (work, ('prepare_packet_patch.py', 'packet-policy.patch', 'partition-proof.json',
            'rows_adapter.cpp', 'rows_bridge.h', 'adapter-compile-result.json', 'compile/adapter/owner.json')),
    (previous, ('provider/source-contract.json',)),
]:
    for name in names:
        path = base / name
        inputs[str(path.relative_to(root))] = pin(path.read_bytes())
assert all(pin((root / name).read_bytes()) == value for name, value in inputs.items())
contract = {
    'scope': 'N packs masked4..7row SIMD tails with full8packets;1..3tails keep independent original F scalar1tasks. Existing CPU/v3 guard retains scalar1allrows on other CPUs. Closed N=8q+r partition proves exact row coverage and unchanged original64*N query counts. Seven D/F/L native objects retain exact paths/bytes; only I adapter3group predicates change under its exact strict g++14argv. Positive B/F guards, L arithmetic/all9380calls/all5programs, separate zero-header no-op, M dead-SSA cleanup, precision/fields/trials/limits/governors/eligibility/ordinaryfallback remain unchanged. Remainder histograms motivate this utilization correction; no full-frame cause/speedup/global GPU attribution or numerical/performance/scientific/release claim.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()), 'source_inputs': inputs,
    'new_native_compilations': 1, 'new_ISPC_compilations': 0,
    'native_objects_reused': 7, 'adapter_recompiled': True,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'medium_tail_partition_proof': pin((work / 'partition-proof.json').read_bytes()),
}
provider.mkdir(exist_ok=True)
(provider / 'provider.patch').write_text(patch)
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
assert before == {str(path.relative_to(source_root)): pin(path.read_bytes())
                  for path in source_root.rglob('*') if path.is_file()}
print(json.dumps({'patch_written': True, 'changed_provider_files': 5,
                  'native_objects_reused': 7, 'new_adapter_compilations': 1,
                  'new_ISPC_compilations': 0, 'shared_source_modified': False}))
