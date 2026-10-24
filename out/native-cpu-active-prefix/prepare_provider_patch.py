"""Emit the exact three-file provider patch for root to review/apply.

This script never changes the shared Mesa source or launches a build.
"""
from pathlib import Path
import difflib
import hashlib
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
basis = root / 'out/native-cpu-slang-bridge'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(records):
    return pin(json.dumps(records, sort_keys=True, separators=(',', ':')).encode())


assert not (provider / 'provider.patch').exists()
compiled = json.loads((work / 'batch-compile-result.json').read_text())
assert len(compiled['objects']) == 3 and all(item['exit'] == 0 for item in compiled['objects'])
owner = json.loads((work / 'build/compile/owner.json').read_text())
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
original = json.loads((basis / 'provider/source-contract.json').read_text())
before = {str(path.relative_to(source_root)): pin(path.read_bytes())
          for path in source_root.rglob('*') if path.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
changes = {
    'src/gallium/targets/lavapipe/meson.build':
        [(str((root / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch' / name).resolve()),
          str((root / 'bin/linux-gcc/native-cpu-active-prefix/typed-batch' / name).resolve()))
         for name in ('transport.o', 'endpoint.o', 'dopri_phase.o')],
    'src/gallium/frontends/lavapipe/lvp_device.c': [('SIR-CPU-25.2.8A', 'SIR-CPU-25.2.8B')],
    'src/gallium/drivers/llvmpipe/lp_screen.c':
        [('Sirius-native-cpu-slang-bridge-v1', 'Sirius-native-cpu-active-prefix-v2')],
}
patch = ''
modified = {}
after = dict(before)
for name, substitutions in changes.items():
    old = (source_root / name).read_text()
    new = old
    for a, b in substitutions:
        assert new.count(a) == 1, (name, a)
        new = new.replace(a, b)
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                         fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(new.encode())
    modified[name] = {'basis': before[name], 'candidate': after[name]}
(provider / 'provider.patch').write_text(patch)
objects = [item['object'] for item in compiled['objects']]
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {k: item[k] for k in ('bytes', 'sha256')}
contract = {'scope': 'Substantive activeprefix bridge repair; original shaderCPP/HPP/SPV/programs and host unchanged. Mesa patch changes only3objectpaths and uniqueprivateBcache/UUID; no arithmetic or original JIT change.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()),
    'source_inputs': {'out/native-cpu-active-prefix/native_packet.hpp': pin((work / 'native_packet.hpp').read_bytes())},
    'limits': 'Configuration/relink/loading/newnativeexecution/C8/ABBA/all8/frame remain unproved.'}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'root_must_review_and_apply': str(provider / 'provider.patch')}))
