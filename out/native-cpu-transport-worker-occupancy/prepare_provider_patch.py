"""Bind the unapplied Q policy to exact M source and all eight reused objects."""
from pathlib import Path
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
objects = original['native_objects']
assert len(objects) == 8
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {key: item[key] for key in ('bytes', 'sha256')}
namespace = runpy.run_path(str(work / 'prepare_patch.py'))
proof = json.loads((work / 'source-proof.json').read_text())
assert not proof['patch_applied'] and proof['new_native_compilations'] == 0
patch = (work / 'provider.patch').read_bytes()
assert pin(patch) == proof['patch']
assert set(namespace['new_text']) == set(proof['changed_files']) and len(namespace['new_text']) == 3
after = dict(before)
for name, text in namespace['new_text'].items():
    assert before[name] == proof['changed_files'][name]['basis']
    after[name] = pin(text.encode())
    assert after[name] == proof['changed_files'][name]['candidate']
inputs = dict(original['source_inputs'])
source_prefix = str(source_root.relative_to(root)) + '/'
for name in proof['changed_files']:
    inputs.pop(source_prefix + name, None)
for name, expected in proof['source_inputs'].items():
    if name.removeprefix(source_prefix) in proof['changed_files']:
        continue
    assert pin((root / name).read_bytes()) == expected, name
    inputs[name] = expected
for name in ('prepare_patch.py', 'provider.patch', 'source-proof.json'):
    path = work / name
    inputs[str(path.relative_to(root))] = pin(path.read_bytes())
for name in ('provider/source-contract.json', 'provider/linked-native-proof.json',
             'control-8-acceptance.json', 'control-8/native-observation.json'):
    path = previous / name
    inputs[str(path.relative_to(root))] = pin(path.read_bytes())
assert all(pin((root / name).read_bytes()) == value for name, value in inputs.items())
contract = {
    'scope': 'Q uses the existing immutable pool worker count only to choose scalar one-row Transport tasks when N fits a nonzero pool. Other positive native stages, W0/larger Transport, zero-active no-op, full B/F guards, ISA policy, DCE, FP/layout/math/counts/governors and all ordinary fallback behavior remain M. All8 exact M objects reused directly; P matched host products reused separately. No latency crossover, numerical/benefit/frame/scientific/release result is claimed.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': proof['changed_files'], 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch), 'source_inputs': inputs,
    'new_native_compilations': 0, 'adapter_recompiled': False,
    'guard_shared_byte_authority': 'out/native-cpu-active-prefix/native_packet.hpp',
    'transport_worker_occupancy_source_proof': pin((work / 'source-proof.json').read_bytes()),
}
provider.mkdir(exist_ok=True)
(provider / 'provider.patch').write_bytes(patch)
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
reuse = {'scope': 'Truthful exact M object reuse receipt; no new compiler result.',
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
