"""Emit the native scheduler correction; reuse all three exact strict objects.

No shared-source mutation, compiler, provider load, or workload is performed.
"""
from pathlib import Path
import difflib
import hashlib
import json

root = Path.cwd()
work = Path(__file__).resolve().parent
provider = work / 'provider'
basis = root / 'out/native-cpu-static-program'
source_root = root / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def canonical(records):
    return pin(json.dumps(records, sort_keys=True, separators=(',', ':')).encode())


assert not (provider / 'provider.patch').exists()
compiled = json.loads((basis / 'batch-compile-result.json').read_text())
owner = json.loads((basis / 'build/compile/owner.json').read_text())
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert compiled['all_prepared_source_bytes_unchanged']
original = json.loads((basis / 'provider/source-contract.json').read_text())
before = {str(path.relative_to(source_root)): pin(path.read_bytes())
          for path in source_root.rglob('*') if path.is_file()}
assert canonical(before) == original['full_candidate_source_identity']
changes = {
    'src/gallium/drivers/llvmpipe/lp_state_cs.c': [
        ('   bool intergroup = llvmpipe->cs->intergroup_row_owner_eligible &&\n',
         '   /* Native scalar groups need one original workgroup per pool task.\n'
         '    * Keep the existing ordinary-ID variant for every packet refusal.\n'
         '    */\n'
         '   bool intergroup = llvmpipe->cs->intergroup_row_owner_eligible &&\n'
         '      !llvmpipe->cs->base.ir.nir->info.cs.llvmpipe_native_cpu_stage &&\n'),
        ('llvmpipe->pipeline_statistics.cs_invocations += num_tasks * info->block[0] * info->block[1] * info->block[2];',
         'llvmpipe->pipeline_statistics.cs_invocations += (uint64_t)num_tasks * info->block[0] * info->block[1] * info->block[2];'),
    ],
    'src/gallium/frontends/lavapipe/lvp_device.c': [('SIR-CPU-25.2.8D', 'SIR-CPU-25.2.8E')],
    'src/gallium/drivers/llvmpipe/lp_screen.c': [
        ('Sirius-native-cpu-static-program-v4', 'Sirius-native-cpu-workgroup-scheduler-v5')],
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
    patch += ''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), n=0,
                                         fromfile='a/' + name, tofile='b/' + name))
    after[name] = pin(new.encode())
    modified[name] = {'basis': before[name], 'candidate': after[name]}
objects = [item['object'] for item in compiled['objects']]
assert len(objects) == 3 and objects == original['native_objects']
for item in objects:
    assert pin((root / item['path']).read_bytes()) == {k: item[k] for k in ('bytes', 'sha256')}
for name, identity in compiled['source_inputs'].items():
    assert pin((root / name).read_bytes()) == identity, name
(provider / 'provider.patch').write_text(patch)
contract = {
    'scope': 'Universal ordinary one-workgroup task schedule for all three exact native stages. Original grid/block counts,64logicalcalls per nativeGroup, all arithmetic/inputs/outputs/guards/pool/FP/wait retained. Eager ordinary JIT key0 is bound before queueing; full packet refusal still uses it. Query multiplication widened before multiplying, preserving counts above32bits. No shape/row/fixture thresholds or worker/governor changes.',
    'basis_source_identity': original['full_candidate_source_identity'],
    'modified_files': modified, 'full_candidate_source_identity': canonical(after),
    'native_objects': objects, 'patch': pin(patch.encode()),
    'source_inputs': compiled['source_inputs'],
    'reused_compile_owner': pin((basis / 'build/compile/owner.json').read_bytes()),
    'reused_objects_not_recompiled': True,
    'limits': 'Source preparation only. Provider build/code association/C8/ABBA/all8/frame/science/release remain unproved.',
}
(provider / 'source-contract.json').write_text(json.dumps(contract, indent=2) + '\n')
print(json.dumps({'patch_written': True, 'shared_source_modified': False,
                  'strict_objects_reused': 3, 'root_must_review_and_apply': str(provider / 'provider.patch')}))
