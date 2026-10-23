"""Authenticate the additive policy against the sealed prior 11,100-file source."""
from pathlib import Path
import hashlib, json, shutil, subprocess, tempfile, zipfile

R = Path.cwd(); W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
def doc(p): return json.loads(p.read_text())
def data_pin(b): return {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
def pin(p): return data_pin(p.read_bytes())
def identity(d): return {k: d[k] for k in ['bytes', 'sha256']}
c = doc(W / 'source-contract.json')
assert pin(W / 'machine-codegen.patch') == identity(c['machine_codegen_patch'])
for name in ['IEEE32_prerequisite', 'prior_intergroup_prototype', 'prior_intergroup_contract', 'prior_ROOT_admission']:
    item = c['source_basis'][name]; assert pin(R / item['path']) == identity(item), name
D = R / 'attestations/qualification/263abe74/lavapipe-intergroup-row-owner-finite-outcome'
archive = D / 'diagnostics-and-trial.zip'
assert pin(archive) == doc(D / 'MANIFEST.json')['archives'][archive.name]
member = 'out/lavapipe-intergroup-row-owner/build/inputs-before.json'
with zipfile.ZipFile(archive) as z: before = z.read(member)
prefix = str(S.relative_to(R)) + '/'
prior = {k[len(prefix):]: identity(v) for k, v in json.loads(before).items() if k.startswith(prefix)}
current = {str(p.relative_to(S)): pin(p) for p in S.rglob('*') if p.is_file()}
assert len(prior) == len(current) == 11100 and prior.keys() == current.keys()
modified = c['modified_files']; assert len(modified) == 8
for name, item in modified.items():
    assert prior[name] == identity(item['basis']), name
    assert current[name] == identity(item['candidate']), name
changed = {name for name in current if current[name] != prior[name]}
assert changed == set(modified), changed
with tempfile.TemporaryDirectory(prefix='source-review-', dir=W) as temporary:
    temporary = Path(temporary)
    for name in modified:
        p = temporary / name; p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(S / name, p)
    subprocess.run(['patch', '--reverse', '--batch', '--fuzz=0', '-p1', '-d', str(temporary),
                    '-i', str(W / 'machine-codegen.patch')], check=True, capture_output=True)
    assert not list(temporary.rglob('*.rej'))
    for name in modified: assert pin(temporary / name) == prior[name], name
soa = 'src/gallium/auxiliary/gallivm/lp_bld_nir_soa.c'
assert current[soa] == {'bytes': 246649, 'sha256': '734425da806076123bdd8fb297a361566db0cd52a70cad9e934a45067b6fbb8e'}
header = (Path('/usr/lib/llvm-21/include/llvm/Config/llvm-config.h')).read_text()
assert '#define LLVM_VERSION_MAJOR 21' in header and '#define LLVM_VERSION_PATCH 8' in header
commands = doc(R / 'bin/linux-gcc/lavapipe-acqrel-ordering/compile_commands.json')
assert commands and all('-DGALLIVM_USE_ORCJIT=0' in row['command'] for row in commands)
for leaf in ['lp_bld_init.c', 'lp_state_cs.c', 'lvp_device.c', 'lp_screen.c', 'lvp_pipeline.c']:
    rows = [row for row in commands if Path(row['file']).name == leaf]
    assert len(rows) == 1 and '-DLP_VULKAN_IEEE32=1' in rows[0]['command'], leaf
cpp = [row for row in commands if Path(row['file']).name == 'lp_bld_misc.cpp']
assert len(cpp) == 1  # Bridge itself is MCJIT/LLVM21 guarded; IEEE32 gate is at its C caller.
host = c['actual_intended_host']; record = host['matched_result']
assert pin(R / record['path']) == identity(record)
matched = doc(R / record['path']); exe = R / matched['executable']
assert pin(exe)['sha256'] == host['executable_sha256'] == matched['executable_sha256']
for kind, record in host['modules'].items():
    module = R / record['module']; assert pin(module) == identity(record)
    assert exe.read_bytes().count(module.read_bytes()) == 1, kind
pipeline = (S / 'src/gallium/frontends/lavapipe/lvp_pipeline.c').read_text()
assert '0xc0,0xf3,0x8b,0x8b' not in pipeline
assert '0xe4,0xf3,0x2b,0xfe' in pipeline and '0xe2,0xbe,0xef,0x14' in pipeline
assert len(c['cache']['Vulkan_UUID_exact16bytes'].encode()) == 16
result = {'passed': True, 'scope': 'Independent ROOT source admission only; machine execution, identical optimized IR, allocator selection and all original gates remain unadmitted.',
    'machine_codegen_patch': pin(W / 'machine-codegen.patch'), 'source_contract': pin(W / 'source-contract.json'),
    'source_files': {name: current[name] for name in sorted(modified)},
    'full_current_source_identity': data_pin(json.dumps(current, sort_keys=True, separators=(',', ':')).encode()),
    'basis': {'preserved_at_commit': '206ca0d4', 'archive': str(archive.relative_to(R)), **pin(archive),
              'member': member, 'member_identity': data_pin(before), 'source_files': len(prior)},
    'changed_files': sorted(changed), 'exact_reverse_recovers_all_eight_prior_file_pins': True,
    'unchanged_remaining_source_files': len(prior) - len(changed),
    'reviewed_contracts': [
        'Only exact E/DP source hashes admit intergroup key1; both Transport hashes use ordinary key0. LocalSize64, original groups/fields/trials and modes remain.',
        'Per-gallivm MCJIT bool is calloc-zero and set only by admitted key1 with IEEE32, MCJIT0 and LLVM21. Header dependency closure must rebuild all struct consumers.',
        'Original engine creation, IR pass call/arguments, CPU/features/FP options/attributes and IEEE32 repair retained. Hook at common cached/uncached join follows original IR passes and precedes emission/publication.',
        'Engine-owned target mutation is mutex-protected, requires Default, sets only CodeGenOptLevel::None; unsupported/unexpected prerequisites refuse before emission. No shared context/global flag/function attribute mutation.',
        'New private UUID/disk namespace plus existing deterministic variant key isolate old/default/None objects. ORC and other LLVM majors are unchanged and unclaimed.',
        'Partial12s stack supports examining codegen but does not establish the300s timeout cause. No speed or compilation completion promise; actual new code and unchanged32/22/all8/frame gates remain required.'
    ], 'restored_original_barrier_lowering': current[soa],
    'original_IEEE32_prerequisite_retained_by_exact_basis_composition': True,
    'new_Sirius_or_shader_compilations': 0, 'runtime_or_benefit_admitted': False}
target = W / 'review-admission.json'; assert not target.exists()
target.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'source_review_passed': True, 'sealed_source_files': len(current), 'additive_changed_files': len(changed),
                  'reused_exact_host_ELF': host['executable_sha256'], 'runtime_or_benefit_admitted': False}))
