"""Preserve the single completed IEEE32 trial without executing any workload."""
from pathlib import Path
import hashlib, io, json, os, subprocess, tarfile, zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
BUNDLE = R / 'attestations/qualification/cb4d62a6/lavapipe-ieee32-rejected'
SOURCE = 'cb4d62a610c199744a7a70d91a1d18932977196f'
ROOTS = [R / 'out' / name for name in [
    'lavapipe-denorm-review', 'lavapipe-native-fp32',
    'lavapipe-fpstate-contract', 'lavapipe-ieee32',
    'lavapipe-ieee32-validation', 'lavapipe-ieee32-ordinary',
    'lavapipe-ieee32-runtime', 'lavapipe-ieee32-readers',
    'lavapipe-ieee32-queue']]
UPSTREAM = R / 'out/lavapipe-ieee32/upstream'

def write(p, value):
    p.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

def identity(p):
    if p.is_symlink():
        target = os.readlink(p)
        return {'type': 'symlink', 'target': target,
                'sha256': hashlib.sha256(target.encode()).hexdigest()}
    with p.open('rb') as f:
        digest = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'type': 'file', 'bytes': p.stat().st_size,
            'sha256': digest, 'mode': p.stat().st_mode & 0o777}

def birth_alive(pid, ticks):
    p = Path('/proc') / str(pid) / 'stat'
    try:
        fields = p.read_text().rsplit(')', 1)[1].split()
        return fields[19] == str(ticks)
    except FileNotFoundError:
        return False

assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not BUNDLE.exists()
BUNDLE.mkdir(parents=True)
owners = {}
for root in ROOTS:
    for p in root.rglob('owner.json'):
        if p.is_relative_to(UPSTREAM):
            continue
        owner = json.loads(p.read_text())
        assert owner['status'] == 'terminal'
        assert owner['remaining_owned_processes'] == [] and not owner['cleanup_errors']
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        for row in owner['observed_births']:
            assert row['absent']
            assert not birth_alive(row['pid'], row['saved_start_ticks'])
        owners[str(p.relative_to(R))] = {
            'identity': identity(p), 'passed': owner['passed'],
            'child_exit': owner['child_exit'], 'stop_reason': owner['stop_reason'],
            'observed_births': owner['observed_births'],
            'wall_seconds': owner['wall_seconds'],
            'peak_sampled_rss_kib': owner['peak_sampled_rss_kib']}
assert len(owners) == 16
write(BUNDLE / 'ownership-closure.json', {
    'all16_terminal_and_reaped': True, 'owners': owners,
    'failed_owners_preserved': [p for p, x in owners.items() if not x['passed']]})

# Authenticate all actual source files against the already pushed original tar
# and the five changed files. No re-extraction, provider build or runtime retry.
admission = json.loads((R / 'out/lavapipe-ieee32/source-admission.json').read_text())
tar_ref = admission['upstream_preservation_revision'] + ':' + admission['official_tar_member']
blob = subprocess.check_output(['git', 'show', tar_ref])
assert hashlib.sha256(blob).hexdigest() == admission['official_tar_sha256']
source_root = UPSTREAM / 'mesa-25.2.8'
source_files = {}
with tarfile.open(fileobj=io.BytesIO(blob), mode='r:xz') as archive:
    for member in archive.getmembers():
        if member.isdir():
            continue
        rel = Path(member.name).relative_to('mesa-25.2.8')
        p = source_root / rel
        if member.issym():
            assert p.is_symlink() and os.readlink(p) == member.linkname
        else:
            assert member.isfile()
            expected = ((R / 'out/lavapipe-ieee32/changed' / rel).read_bytes()
                        if str(rel) in admission['changed_files']
                        else archive.extractfile(member).read())
            assert p.read_bytes() == expected, rel
        source_files[str(p.relative_to(R))] = identity(p)
actual_files = {str(p.relative_to(R)) for p in source_root.rglob('*')
                if p.is_symlink() or p.is_file()}
assert actual_files == set(source_files)
write(BUNDLE / 'upstream-reconstruction.json', {
    'original_tar_git_object': tar_ref,
    'original_tar_sha256': admission['official_tar_sha256'],
    'candidate_patch': 'out/lavapipe-ieee32/candidate.patch',
    'all_actual_source_paths_reconstruct_exactly': True,
    'files': source_files})

def pack(destination, files):
    manifest = {}
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED,
                         compresslevel=6, allowZip64=True) as archive:
        for p in sorted(set(files)):
            assert p.is_file() and not p.is_symlink(), p
            rel = str(p.relative_to(R))
            manifest[rel] = identity(p)
            archive.write(p, rel)
    # Verify every original byte through the actual newly created container.
    with zipfile.ZipFile(destination) as archive:
        assert set(archive.namelist()) == set(manifest)
        assert archive.testzip() is None
        for name, record in manifest.items():
            content = archive.read(name)
            assert len(content) == record['bytes']
            assert hashlib.sha256(content).hexdigest() == record['sha256']
    return {'archive': identity(destination), 'files': manifest}

diagnostics = [p for root in ROOTS for p in root.rglob('*')
               if p.is_file() and not p.is_relative_to(UPSTREAM)]
diagnostics.append(Path(__file__).resolve())
products = [p for name in ['lavapipe-native-fp32', 'lavapipe-fpstate-contract',
                           'lavapipe-ieee32', 'lavapipe-ieee32-controls',
                           'lavapipe-ieee32-staged']
            for p in (R / 'bin/linux-gcc' / name).rglob('*') if p.is_file()]
actions = json.loads((R / 'out/lavapipe-ieee32-runtime/actions.json').read_text())
products += [R / p for a in actions.values() for p in a['frozen_inputs']
             if p.startswith('bin/linux-gcc/src/')]
containers = {
    'diagnostics-and-trial.zip': pack(BUNDLE / 'diagnostics-and-trial.zip', diagnostics),
    'providers-and-products.zip': pack(BUNDLE / 'providers-and-products.zip', products)}

benefit = json.loads((R / 'out/lavapipe-ieee32-readers/benefit.json').read_text())
decimal = json.loads((R / 'out/lavapipe-ieee32-readers/decimal-agreement.json').read_text())
assert benefit['audit_passed'] and not benefit['benefit_passed']
assert benefit['checks_total'] == 22 and benefit['checks_passed'] == 11
assert decimal['agrees_on_all22'] and not decimal['all22_pass']
assert not any((R / 'out/lavapipe-ieee32-runtime' / f'control-{i}').exists()
               for i in range(1, 8))
checks = benefit['checks'] + benefit['accounting_checks']
summary = {
    'source': SOURCE, 'public_source': '0109ab462edade524667d20303507aa09adcbaf6',
    'decision': 'Rejected for performance adoption; real provider repair remains a private finite prototype.',
    'normal_RTE_OpFAdd_counterexample': {
        'row': 1508, 'a': '0x00800000', 'b': '0x3f000000',
        'expected': '0x3f000000', 'actual_baseline_poisoned_caller': '0x3f000001'},
    'actual_provider_repair': 'Scoped nearest-even, gradual-underflow execution and full caller/worker MXCSR restoration; five Mesa source files in private Vulkan-only x86-64 build.',
    'primitive_acceptance': 'Both original candidate default/caller-poison cases: 174268 records, all522804 words exact, observed callback states restored.',
    'ordinary_native_Transport': {'actual_called_function_bytes': 336181,
        'packed_VDIVPS': 9, 'packed_VSQRTPS': 1, 'all9_original_memory_fences': True},
    'C8_pass': True, 'sole_ABBA22': {
        'audit_pass': True, 'all32_complete_readback_digests_cross_arm_equal': True,
        'total': 22, 'passed': 11,
        'whole_call_passed': sum(c['pass'] for c in checks if c['scope'] == 'complete_call_ms'),
        'stages_passed': sum(c['pass'] for c in checks if 'inclusive Dispatch' in c['scope']),
        'provisioning_API_passed': sum(c['pass'] for c in benefit['accounting_checks']),
        'decimal_agrees_on_all22': True, 'checks': checks},
    'remaining7_executed': False, 'cold_frame_executed': False,
    'unchanged_retries': 0, 'full_scientific_or_Mandatory_or_release_qualification': False,
    'issue': 'https://github.com/sunyaisblank/Sirius/issues/106',
    'next_work_basis': 'Use actual stage costs and source ownership to eliminate duplicate native work before another candidate; preserve original arithmetic, workgroups, all gates and full call accounting.'}
write(BUNDLE / 'findings.json', summary)
write(BUNDLE / 'MANIFEST.json', {
    'schema': 'sirius-private-ieee32-single-rejected-trial-v1',
    'source': SOURCE, 'adopted': False, 'full_qualification': False,
    'containers': containers,
    'other_files': {p.name: identity(p) for p in BUNDLE.iterdir()
                    if p.is_file() and p.suffix != '.zip'}})
print(json.dumps({'bundle': str(BUNDLE.relative_to(R)),
                  'owners': len(owners),
                  'containers': {k: v['archive'] for k, v in containers.items()},
                  'source_paths': len(source_files)}, indent=2))
