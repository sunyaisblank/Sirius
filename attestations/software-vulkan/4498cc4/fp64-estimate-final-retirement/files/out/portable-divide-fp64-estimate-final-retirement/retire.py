"""Retire only completed, pushed and byte-matched private experiment roots."""
from pathlib import Path
import hashlib, json, shutil, subprocess

R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
query = lambda *args: subprocess.check_output(list(args), cwd=R, text=True).strip()
assert query('git', 'rev-parse', 'HEAD') == '4498cc41168ad29fc4c111819a4cb7ccaa54a8bf'
assert not query('git', 'status', '--porcelain=v1')
assert not (W / 'retirement.json').exists()
groups = [
    ('fp64-estimate-original-arithmetic', 'portable-divide-fp64-estimate-arithmetic', True),
    ('fp64-estimate-controls-profile-reader-stop', 'portable-divide-fp64-estimate-controls', True),
    ('fp64-estimate-profile-recovery', 'portable-divide-fp64-estimate-profile-recovery', False),
    ('fp64-estimate-original-controls', 'portable-divide-fp64-estimate-controls-continuation', True),
    ('fp64-estimate-original-benefit', 'portable-divide-fp64-estimate-benefit', True),
]

def content(path):
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}

def inventory(path):
    assert path.is_dir() and not path.is_symlink()
    entries = sorted(path.rglob('*'))
    assert all(not entry.is_symlink() for entry in entries)
    return {str(entry.relative_to(path)): content(entry) for entry in entries if entry.is_file()}

refs = ['refs/tags/evidence/issue-44-' + name for name, _, _ in groups]
remote = dict((row.split()[1], row.split()[0]) for row in query('git', 'ls-remote', 'origin', *refs).splitlines())
rows, checked_births = [], set()
for bundle_name, private_name, has_products in groups:
    bundle = R / 'attestations/software-vulkan/4bcd987' / bundle_name
    tag = 'refs/tags/evidence/issue-44-' + bundle_name
    commit = query('git', 'rev-parse', tag)
    assert remote[tag] == commit
    manifest = json.loads((bundle / 'manifest.json').read_text())
    if bundle_name.endswith('original-benefit'):
        disposition = json.loads((bundle / 'files/out' / private_name / 'disposition.json').read_text())
        assert disposition['decision'] == 'REJECT_GUARDED_FP64_ESTIMATE_COST'
        assert disposition['original_benefit_trial_completed'] and not disposition['production_adopted']
    paths = ['out/' + private_name]
    if has_products:
        paths.append('bin/linux-gcc/' + private_name)
    if bundle_name == 'fp64-estimate-original-controls':
        paths.append('out/original-cold-kerr-source-review-continuation')
    for relative in paths:
        live = R / relative
        saved = bundle / 'files' / relative
        assert not query('git', 'ls-files', '--', relative)
        before = inventory(live)
        assert before == inventory(saved), relative
        for name, pin in before.items():
            member = 'files/' + relative + '/' + name
            assert manifest['files'][member] == pin
            expected_blob = query('git', 'rev-parse', commit + ':' + str((bundle / member).relative_to(R)))
            assert query('git', 'hash-object', '--', str(live / name)) == expected_blob
        owners = list(live.rglob('owner.json'))
        for owner_path in owners:
            owner = json.loads(owner_path.read_text())
            assert owner['status'] == 'terminal' and owner['cleanup_errors'] == []
            assert owner['remaining_owned_processes'] == [] and owner['observed_births']
            for birth in owner['observed_births']:
                assert birth['absent']
                stat = Path('/proc') / str(birth['pid']) / 'stat'
                try:
                    current = int(stat.read_text().split(') ', 1)[1].split()[19])
                except FileNotFoundError:
                    current = None
                assert current != birth['saved_start_ticks']
                checked_births.add((birth['pid'], birth['saved_start_ticks']))
        rows.append({'root': relative, 'bundle': str(bundle.relative_to(R)), 'commit': commit,
                     'manifest': content(bundle / 'manifest.json'), 'files': len(before),
                     'file_map_sha256': hashlib.sha256(json.dumps(before, sort_keys=True).encode()).hexdigest(),
                     'terminal_owner_records': len(owners), 'inventory': before})

assert len(rows) == 10 and len({row['root'] for row in rows}) == 10
removed = []
try:
    for row in rows:
        live = R / row['root']
        assert inventory(live) == row['inventory']
        shutil.rmtree(live)
        assert not live.exists()
        removed.append(row['root'])
finally:
    report = {'schema': 'sirius-fp64-estimate-completed-private-retirement-v1',
              'completed': len(removed) == 10, 'removed': removed,
              'roots': [{key: value for key, value in row.items() if key != 'inventory'} for row in rows],
              'fresh_recorded_birth_absences': len(checked_births),
              'source_revision': query('git', 'rev-parse', 'HEAD'),
              'tracked_status_after': query('git', 'status', '--porcelain=v1'),
              'scope': 'Ten completed private roots matched pushed immutable Git/manifest bytes before removal. Original failed reader and fixed rejected performance verdict remain preserved. No full qualification or systemwide consumer-absence claim.'}
    (W / 'retirement.json').write_text(json.dumps(report, sort_keys=True, indent=2) + '\n')
assert report['completed'] and report['tracked_status_after'] == ''
print(json.dumps({'completed': True, 'removed_roots': len(removed),
                  'fresh_recorded_birth_absences': len(checked_births)}))
