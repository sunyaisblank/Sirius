"""Bind ROOT's source review to the preserved basis and actual shader products."""
from pathlib import Path
import hashlib, io, json, shutil, subprocess, tarfile, tempfile, zipfile

R = Path.cwd(); W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
def doc(p): return json.loads(p.read_text())
def data_pin(b): return {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
def pin(p): return data_pin(p.read_bytes())
def identity(d): return {k: d[k] for k in ('bytes', 'sha256')}
c = doc(W / 'source-contract.json')
assert pin(W / 'prototype.patch') == identity(c['prototype_patch'])
assert pin(R / c['root_owned_reused_products']['path']) == identity(c['root_owned_reused_products'])
assert pin(R / c['basis']['IEEE32_prerequisite']['path']) == identity(c['basis']['IEEE32_prerequisite'])
archive = R / 'attestations/qualification/1c95c98f/lavapipe-acqrel-ordering-finite-outcome/diagnostics-and-trial.zip'
assert pin(archive) == {'bytes': 7634909, 'sha256': 'a29355e06f4d6d6e16971ecd1868d855d53e434088669fe4631ddf5f3ee2baf9'}
member = 'out/lavapipe-acqrel-ordering/build/inputs-before.json'
with zipfile.ZipFile(archive) as z: prior_bytes = z.read(member)
prefix = str(S.relative_to(R)) + '/'
prior = {k[len(prefix):]: identity(v) for k, v in json.loads(prior_bytes).items() if k.startswith(prefix)}
current = {str(p.relative_to(S)): pin(p) for p in S.rglob('*') if p.is_file()}
assert len(current) == len(prior) == 11100
assert current.keys() == prior.keys()
modified = c['modified_files']
assert len(modified) == 6
for name, d in modified.items():
    assert prior[name] == identity(d['basis']), name
    assert current[name] == identity(d['candidate']), name
with tempfile.TemporaryDirectory(prefix='source-review-', dir=W) as temporary:
    temporary = Path(temporary)
    for name in modified:
        p = temporary / name; p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(S / name, p)
    reverse = subprocess.run(['patch', '--reverse', '--batch', '--fuzz=0', '-p1', '-d', str(temporary),
                              '-i', str(W / 'prototype.patch')], check=True, capture_output=True, text=True)
    assert not list(temporary.rglob('*.rej'))
    for name in modified: assert pin(temporary / name) == prior[name], name
restored = c['restored_original_barrier_lowering']
restored_name = str(Path(restored['path']).relative_to(S.relative_to(R)))
release = c['basis']['release_archive']
release_bytes = subprocess.check_output(['git', 'show', release['git']])
assert data_pin(release_bytes) == identity(release)
with tarfile.open(fileobj=io.BytesIO(release_bytes), mode='r:xz') as tar:
    original_soa = tar.extractfile('mesa-25.2.8/' + restored_name).read()
assert data_pin(original_soa) == identity(restored) == current[restored_name]
actual_changes = {name for name in current if current[name] != prior[name]}
assert actual_changes == set(modified) | {restored_name}, actual_changes
for name, d in c['basis']['IEEE32_files'].items():
    assert prior[name] == identity(d['basis'])
    assert current[name] == identity(d['current'])
reused = doc(R / c['root_owned_reused_products']['path'])
for name, d in reused['restored_actual_products'].items():
    p = R / name
    assert pin(p) == identity(d) and (p.stat().st_mode & 0o777) == d['mode'], name
for name, d in reused['reused_original_static_and_code_proof'].items(): assert pin(R / name) == identity(d), name
for name, d in c['eligible_modules'].items():
    p = R / name; assert pin(p) == identity(d)
    assert hashlib.sha1(p.read_bytes()).hexdigest() == d['sha1_for_original_Vulkan_cache_identity']
    lines = p.with_suffix('.spvasm').read_text().splitlines()
    entries = [line.strip() for line in lines if 'OpEntryPoint ' in line]
    assert len(entries) == 1 and 'OpEntryPoint GLCompute %ComputeMain "main"' in entries[0]
    assert 'OpExecutionMode %ComputeMain LocalSize 64 1 1' in '\n'.join(lines)
assert '.pName = "main"' in (R / 'src/sirius/backend/vulkan/vulkan_device.cpp').read_text()
assert 'strcmp(sinfo->pName, "main")' in (S / 'src/gallium/frontends/lavapipe/lvp_pipeline.c').read_text()
result = {
    'passed': True, 'scope': 'Independent ROOT source admission only; actual generated code and unchanged runtime gates remain required.',
    'prototype_patch': pin(W / 'prototype.patch'), 'source_contract': pin(W / 'source-contract.json'),
    'source_files': {name: current[name] for name in sorted(set(modified) | set(c['basis']['IEEE32_files']) | {restored_name})},
    'preserved_prior_source': {'commit': 'd324d909463615699517b9aa5034bb521ceddbd2',
        'archive': str(archive.relative_to(R)), **pin(archive), 'member': member,
        'member_identity': data_pin(prior_bytes), 'source_files': len(prior)},
    'full_current_source_identity': data_pin(json.dumps(current, sort_keys=True, separators=(',', ':')).encode()),
    'changed_from_prior_trial': sorted(actual_changes),
    'reverse_patch_recovers_exact_prior_six_file_basis': True,
    'original_barrier_source_authenticated_from_release': True,
    'original_IEEE32_prerequisite_retained': True,
    'original_spirv_entrypoint_is_main_not_debug_ComputeMain': True,
    'reviewed_contracts': [
        'Only three exact immutable modules; one function, local0-only observable owner; unsupported source/shape/base/width and binary imports fall back.',
        'Variant-local options make group IDs divergent; singular-function gate prevents scalar helper-call context, and lp_build_opt_nir has no eligible helper-introducing pass.',
        'Eight X-bounded independent group owners; original Y/Z/grid/block/query counts; out-of-grid mask precedes the shader prologue and observable accesses.',
        'Private state, task-pool completion and original scoped IEEE32 save/set/restore retained; failed acquire-release modification removed.',
        'Private binary UUID, disk-cache salt, serialized shader marker and separate variant key prevent stale/ordinary/intergroup cache aliasing.'
    ],
    'prebuild_corrections': ['Vulkan entrypoint main verified against all three immutable OpEntryPoint declarations and the actual Vulkan pipeline factory.'],
    'runtime_or_performance_admitted': False, 'new_builds_or_workloads': 0
}
target = W / 'review-admission.json'; assert not target.exists()
target.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'source_review_passed': True, 'sealed_source_files': len(current),
                  'changed_from_prior_trial': len(actual_changes), 'new_builds_or_workloads': 0}))
