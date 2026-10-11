"""Preserve the rejected trial's final records, then remove its finite staging."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
assert ROOT == Path('/home/astra/.project/Sirius')
assert WORK == ROOT / 'out/native-narrow-fma-trial'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def identity(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def expected(v): return {k: v[k] for k in ('bytes', 'sha256')}
def key(v): return v['bytes'], v['sha256']
def git(*args): return subprocess.check_output(['git', *args], text=True).strip()
source = read(WORK / 'restoration/source.json')
head = source['source_revision']
def source_check():
    assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}')
    assert not git('status', '--porcelain')
    assert git('rev-parse', 'HEAD^{tree}') == git('rev-parse', source['baseline_revision'] + '^{tree}')
    assert git('rev-parse', 'main') == git('rev-parse', 'origin/main') == '5e706073a3fd2f9afcf8efbcde2f094122050f6a'
    assert len([s for s in git('worktree', 'list', '--porcelain').splitlines() if s.startswith('worktree ')]) == 1
source_check()
assert read(WORK / 'native-disposition-reviewed.json')['reject_candidate']
assert read(WORK / 'native-disposition-reviewed.json')['failed_checks'] == 25
assert read(WORK / 'restoration-independent-review.json')['accepted']
for case in ('configure', 'build', 'products', 'bindings', 'admission', 'selection'):
    owner = read(WORK / 'restoration' / case / 'owner.json')
    assert owner['passed'] and owner['child_exit'] == 0 and owner['source'] == {'revision': head, 'status': ''}
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
ci = read(WORK / 'restoration/ci-final.json')
assert ci['headSha'] == head and ci['status'] == 'completed' and ci['conclusion'] == 'success'
jobs = {j['name']: j for j in ci['jobs']}
assert set(jobs) == {'integration-no-render', 'integration-windows-no-render', 'integration-macos-no-render', 'windows-build', 'macos-build', 'linux-sanitizers', 'linux-gate'}
assert all(jobs[n]['conclusion'] == 'success' for n in ('integration-no-render', 'integration-windows-no-render', 'integration-macos-no-render'))
assert all(jobs[n]['conclusion'] == 'skipped' for n in ('windows-build', 'macos-build', 'linux-sanitizers', 'linux-gate'))
clearance = read(WORK / 'cleanup-source-review.json')
assert clearance['accepted'] and clearance['source_revision'] == head
for name, v in clearance['controller_seals'].items(): assert identity(WORK / name) == v
preserved = read(WORK / 'closing-tools-preservation.json')
assert preserved['pushed'] and preserved['source_revision'] == head
assert git('rev-parse', preserved['tag'] + '^{commit}') == preserved['commit']
assert git('ls-remote', 'origin', 'refs/tags/' + preserved['tag']).split()[0] == preserved['commit']
recoveries = {}
verified_recoveries = set()
def verify_recovery(v, r):
    marker = (key(v), json.dumps(r, sort_keys=True))
    if marker in verified_recoveries: return
    if r['kind'] == 'immutable-Git-blob':
        b = subprocess.check_output(['git', 'cat-file', 'blob', r['blob']])
        assert {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()} == expected(v)
    else:
        p = ROOT / r['path']
        assert p.is_file() and not p.is_symlink() and identity(p) == expected(v)
    verified_recoveries.add(marker)
for name, review_name in (('trial-protection.json', 'trial-preservation-reviewed.json'), ('restoration-protection.json', 'restoration-preservation-reviewed.json')):
    receipt = read(WORK / name); review = read(WORK / review_name)
    assert review['accepted'] and review['manifest'] == receipt['manifest']
    archive = ROOT / receipt['archive']
    assert identity(archive / 'manifest.json') == receipt['manifest']
    manifest = read(archive / 'manifest.json')
    for v in manifest['payloads']:
        assert identity(archive / v['path']) == expected(v)
    groups = ('executed_input_ledger', 'diagnostics', 'volatile_product_and_native_stage_files') if name == 'trial-protection.json' else ('executed_input_ledger', 'records')
    for group in groups:
        for v in manifest[group]:
            verify_recovery(v, v['recovery']); recoveries[key(v)] = v['recovery']
for v in preserved['files']:
    assert identity(ROOT / v['origin']) == expected(v)
    r = {'kind': 'immutable-Git-blob', 'revision': preserved['commit'], 'blob': v['blob'], 'path': v['path']}
    assert git('rev-parse', preserved['commit'] + ':' + v['path']) == v['blob']
    verify_recovery(v, r); recoveries[key(v)] = r

plan = read(WORK / 'cleanup-targets.json')
assert plan['source_revision'] == head
allowed_stages = {f'bin/{preset}/native-{rev[:7]}' for preset in ('windows-msvc', 'macos-clang') for rev in (source['baseline_revision'], source['rejected_revision'])}
assert set(plan['remove_stages']) | set(plan['retain_stages']) == allowed_stages
assert not (set(plan['remove_stages']) & set(plan['retain_stages']))
assert all(v['reason'] for v in plan['retained_dependencies'])
assert {v['path'] for v in plan['retained_dependencies']} == set(plan['retain_stages'])
stages = [ROOT / name for name in plan['remove_stages']]
extras = [ROOT / 'bin/linux-gcc/src/sirius/backend/retained' / name for name in ('retained_endpoint_narrow_fma.spv', 'retained_transport_narrow_fma.spv')]
resources = [ROOT / 'bin/linux-gcc/tests/backend/resources', ROOT / 'bin/linux-gcc/src/sirius/app/resources']
def entries(folder):
    assert folder.is_dir() and not folder.is_symlink()
    result = {}
    for p in folder.rglob('*'):
        assert not p.is_symlink() and (p.is_file() or p.is_dir())
        result[str(p.relative_to(folder))] = 'file' if p.is_file() else 'directory'
    return result
def tree_seals(folder):
    entries(folder)
    return {str(p.relative_to(ROOT)): identity(p) for p in sorted(folder.rglob('*')) if p.is_file()}
resource_seals = {str(p.relative_to(ROOT)): tree_seals(p) for p in resources}
assert [len(v) for v in resource_seals.values()] == [23, 32]
retained_stage_seals = {name: tree_seals(ROOT / name) for name in plan['retain_stages']}
assert all(len(v) == 94 for v in retained_stage_seals.values())
product_seals = tree_seals(ROOT / 'bin/linux-gcc/src/sirius/backend/retained')
assert len(product_seals) == 37
assert sum(name.endswith('.spv') for name in product_seals) == 35
assert any(name.endswith('/retained_dense_program_identity.slang') for name in product_seals)
publication = ROOT / 'attestations/review' / head[:7] / 'native-narrow-fma-publication'
cleanup = ROOT / 'attestations/cleanup' / head[:7] / 'native-narrow-fma-finished'
assert not publication.exists() and not cleanup.exists()
deletions, payloads, directory_sets = [], {}, {}
for folder in (WORK, *stages):
    directory_sets[str(folder)] = entries(folder)
    if folder != WORK: assert sum(v == 'file' for v in directory_sets[str(folder)].values()) == 94
    for p in sorted(folder.rglob('*')):
        if not p.is_file(): continue
        v = identity(p); r = recoveries.get(key(v))
        if r is None:
            name = 'payloads/' + v['sha256']; q = publication / name
            q.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(p, q)
            assert identity(q) == v
            payloads[name] = {'path': name, **v}
            r = {'kind': 'local-publication-payload', 'path': str(q.relative_to(ROOT))}
            recoveries[key(v)] = r
        verify_recovery(v, r)
        deletions.append({'origin': str(p.relative_to(ROOT)), **v, 'recovery': r})
for p in extras:
    assert p.is_file() and not p.is_symlink()
    v = identity(p); r = recoveries[key(v)]; verify_recovery(v, r)
    deletions.append({'origin': str(p.relative_to(ROOT)), **v, 'recovery': r})

known = {}
boot_id = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
def birth(pid):
    try: return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError): return None
def absent(pid, ticks):
    assert type(pid) is int and pid > 0 and ticks
    current = birth(pid); assert current != str(ticks), (pid, ticks)
    known[(pid, str(ticks))] = {'pid': pid, 'saved_start_ticks': str(ticks), 'current_start_ticks': current, 'absent': True}
for p in WORK.rglob('owner.json'):
    o = read(p)
    if 'observed_births' not in o: continue
    assert o['boot_id'] == boot_id and not o['cleanup_errors'] and not o['remaining_owned_processes']
    absent(o['controller']['pid'], o['controller']['start_ticks'])
    for v in o['observed_births']: absent(v['pid'], v['saved_start_ticks'])
for p in WORK.glob('*-bridge.json'):
    o = read(p)
    assert not o.get('errors') and not o.get('outer_stop')
    for pid_key, ticks_key in (('linux_owner_pid', 'linux_owner_start_ticks'), ('linux_bridge_child_pid', 'linux_bridge_child_start_ticks'), ('linux_child_pid', 'linux_child_birth')):
        if pid_key in o: absent(o[pid_key], o[ticks_key])
for p in WORK.glob('native-*-sequence.json'):
    o = read(p); assert o['completed']; absent(o['linux_owner_pid'], o['linux_owner_start_ticks'])
    for v in o['owned_calls']:
        assert v['returncode'] == 0 and v['ownership_verified']; absent(v['pid'], v['start_ticks'])
windows = read(WORK / 'fresh-native-consumers.json')
bridge = read(WORK / 'fresh-native-consumers-bridge.json')
accepted = read(WORK / 'fresh-native-consumers-accepted.json')
assert windows['pass'] and not windows['matches'] and len(windows['known_births']) == 45
assert all(v['owned_birth_absent'] for v in windows['known_births'])
assert set(windows['targets']) == {'native-narrow-fma-trial', 'native-b81e061', 'native-21e55ce'}
assert bridge['accepted'] and bridge['returncode'] == 0 and not bridge['errors'] and not bridge['outer_stop'] and not bridge['controller_signals']
assert accepted['pass'] and accepted['terminal_receipt_sha256'] == identity(WORK / 'fresh-native-consumers.json')['sha256']
assert accepted['terminal_bridge_sha256'] == identity(WORK / 'fresh-native-consumers-bridge.json')['sha256']
def fresh_windows():
    age = datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(windows['checked_utc'].replace('Z', '+00:00'))
    assert datetime.timedelta() <= age <= datetime.timedelta(seconds=60), 'Refresh the Windows census immediately before cleanup'
fresh_windows()
targets = ['native-narrow-fma-trial', *[p.name for p in stages], *[str(p.relative_to(ROOT)) for p in extras]]
matches, errors, count = [], [], 0
self_birth = birth(os.getpid())
for p in Path('/proc').glob('[0-9]*'):
    if int(p.name) == os.getpid():
        assert birth(os.getpid()) == self_birth; continue
    count += 1
    for field in ('cmdline', 'environ', 'maps', 'cwd', 'exe'):
        try:
            value = os.readlink(p / field) if field in ('cwd', 'exe') else (p / field).read_bytes().decode(errors='replace')
            if any(target in value for target in targets): matches.append({'pid': int(p.name), 'field': field})
        except (OSError, ValueError): errors.append({'pid': int(p.name), 'field': field})
assert not matches, matches
source_check()
for v in deletions:
    assert identity(ROOT / v['origin']) == expected(v); verify_recovery(v, v['recovery'])
for folder in (WORK, *stages): assert entries(folder) == directory_sets[str(folder)]
assert {str(p.relative_to(ROOT)): tree_seals(p) for p in resources} == resource_seals
publication.mkdir(parents=True, exist_ok=True)
manifest = {'source_revision': head, 'payloads': list(payloads.values()), 'payload_count': len(payloads), 'payload_bytes': sum(v['bytes'] for v in payloads.values()), 'deletion_recovery_ledger': deletions, 'resource_seals': resource_seals, 'retained_stage_seals': retained_stage_seals, 'closing_tag': preserved, 'scope': 'Exact final owner/controller/review/CI/publication/census records and full deletion recovery ledger; prior immutable archives and pushed source tags referenced. Ignored local evidence, not a remote binary backup or scientific qualification.'}
(publication / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
cleanup.mkdir(parents=True)
receipt = {'source_revision': head, 'source_committed_pushed': True, 'publication_archive': str(publication.relative_to(ROOT)), 'publication_manifest': identity(publication / 'manifest.json'), 'removed': False, 'files': len(deletions), 'bytes': sum(v['bytes'] for v in deletions), 'targets': [str(p.relative_to(ROOT)) for p in (WORK, *stages, *extras)], 'retained_dependencies': plan['retained_dependencies'], 'CI': {'run': ci['databaseId'], 'headSha': head, 'status': ci['status'], 'conclusion': ci['conclusion'], 'full_qualification_jobs_skipped': 4}, 'Linux_known_births': list(known.values()), 'Linux_accessible_processes': count, 'Linux_matches': matches, 'Linux_read_errors': errors, 'Linux_self_exclusion': {'pid': os.getpid(), 'start_ticks': self_birth}, 'Windows_scoped_census': windows, 'resources_files_unchanged': [23,32], 'scope': 'Finished rejected trial/restoration scratch, plan-selected reconstructed native stages and two undeclared shader byproducts only. Official native exports, raw/input archives, canonical resources/current original products, tools/build trees retained. Finite accessible scans; inaccessible fields/handles/global/scheduled consumers unverified. No scientific/render/release admission.'}
(cleanup / 'cleanup.json').write_text(json.dumps(receipt, indent=2) + '\n')
fresh_windows()
for p in extras: p.unlink()
for folder in (WORK, *stages): shutil.rmtree(folder)
assert not any(p.exists() for p in (WORK, *stages, *extras))
assert {str(p.relative_to(ROOT)): tree_seals(p) for p in resources} == resource_seals
assert {name: tree_seals(ROOT / name) for name in plan['retain_stages']} == retained_stage_seals
for name, v in product_seals.items():
    if ROOT / name not in extras: assert identity(ROOT / name) == v
assert len(tree_seals(ROOT / 'bin/linux-gcc/src/sirius/backend/retained')) == 35
assert len(list((ROOT / 'bin/linux-gcc/src/sirius/backend/retained').glob('*.spv'))) == 33
source_check()
receipt.update(removed=True, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), out_empty=not any((ROOT / 'out').iterdir()), resources_verified_after=True, original_products_verified_after=True, final_tracking=git('rev-list','--left-right','--count','HEAD...@{u}'))
(cleanup / 'cleanup.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'cleanup': str(cleanup.relative_to(ROOT)), 'files': receipt['files'], 'bytes': receipt['bytes'], 'out_empty': receipt['out_empty'], 'publication_manifest': receipt['publication_manifest']}))
