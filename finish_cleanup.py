"""Recover and remove only the finished #94 scratch directory, without following links."""
from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess

ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
assert WORK == ROOT / 'out/guarded-upmultiply-review'
read = lambda p: json.loads(p.read_text())
def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def expected(v): return {k: v[k] for k in ('bytes', 'sha256')}
def write(path, value): path.write_text(json.dumps(value, indent=2) + '\n')
head = read(WORK / 'restoration/source.json')['source_revision']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert subprocess.check_output(['git', 'rev-parse', 'HEAD']) == subprocess.check_output(['git', 'rev-parse', '@{u}'])
assert not read(WORK / 'software-comparison.json')['pass']
assert read(WORK / 'trial-preservation-reviewed.json')['pass']
assert read(WORK / 'restoration-preservation-reviewed.json')['accepted']
assert read(WORK / 'restoration-source-independent-review.json')['pass']
assert read(WORK / 'cleanup-source-review.json')['pass']
assert read(WORK / 'restoration/build-and-arrays.json')['original40_complete_arrays_and_whole_header_byte_exact95']
assert read(WORK / 'restoration/linux-readonly-payloads.json')['pass_']
ci = read(WORK / 'restoration-ci.json')
assert ci['headSha'] == head and ci['status'] in ('queued', 'in_progress', 'completed')
if ci['status'] == 'completed': assert ci['conclusion'] == 'success'
tools = read(WORK / 'tools-preservation-v3.json')
assert tools['pushed'] and tools['source_revision'] == head
assert subprocess.check_output(['git', 'rev-parse', tools['tag'] + '^{commit}'], text=True).strip() == tools['commit']
assert subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + tools['tag'] + '^{}'], text=True).split()[0] == tools['commit']
for name, value in tools['files'].items():
    data = subprocess.check_output(['git', 'show', tools['tag'] + ':' + name])
    assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected(value)
    assert seal(WORK / value.get('origin', name)) == expected(value)

retained = read(WORK / 'retained-paths.json')
def check_retained():
    for path, value in retained['files'].items():
        p = ROOT / path
        assert p.is_file() and not p.is_symlink() and seal(p) == expected(value)
check_retained()
resources = WORK / 'resources'
target = ROOT / 'bin/linux-gcc/tests/backend/resources'
assert resources.is_symlink() and resources.resolve() == target.resolve()
resource_paths = {str(p.relative_to(ROOT)) for p in target.rglob('*') if p.is_file()}
assert len(resource_paths) == 23 and resource_paths <= set(retained['files'])

recoveries = {}
for name in ('trial-protection.json', 'restoration-protection.json'):
    receipt = read(WORK / name)
    archive = ROOT / receipt['archive']
    assert seal(archive / 'manifest.json') == receipt['manifest']
    manifest = read(archive / 'manifest.json')
    assert {str(p.relative_to(archive)) for p in archive.rglob('*') if p.is_file()} == {'manifest.json', *[v['path'] for v in manifest['payloads']]}
    for value in manifest['payloads']:
        p = archive / value['path']
        assert p.is_file() and not p.is_symlink() and seal(p) == expected(value)
        recoveries.setdefault(value['origin'], []).append({**value, 'recovery': str(p.relative_to(ROOT))})

def entries():
    result = {}
    for directory, dirs, files in os.walk(WORK, followlinks=False):
        for name in [*dirs, *files]:
            p = Path(directory) / name
            relative = str(p.relative_to(WORK))
            if p.is_symlink():
                assert p == resources
                result[relative] = {'kind': 'symlink', 'target': os.readlink(p)}
            elif p.is_dir(): result[relative] = {'kind': 'directory'}
            else:
                assert p.is_file()
                result[relative] = {'kind': 'file', **seal(p)}
    return result
tree = entries()
publication = ROOT / 'attestations/review' / head[:7] / 'rejected-upmultiply-publication'
cleanup = ROOT / 'attestations/cleanup' / head[:7] / 'rejected-upmultiply-finished'
assert not publication.exists() and not cleanup.exists()
payloads, deletions = [], []
for relative, entry in sorted(tree.items()):
    if entry['kind'] != 'file': continue
    path = WORK / relative
    origin = str(path.relative_to(ROOT))
    value = expected(entry)
    recovery = next((v['recovery'] for v in recoveries.get(origin, []) if expected(v) == value), None)
    if recovery is None:
        output = publication / 'diagnostic' / relative
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, output)
        assert seal(output) == value
        payloads.append({'path': str(output.relative_to(publication)), 'origin': origin, **value})
        recovery = str(output.relative_to(ROOT))
    deletions.append({'origin': origin, **value, 'recovery': recovery})
publication.mkdir(parents=True, exist_ok=True)
write(publication / 'manifest.json', {'revision': head, 'payloads': payloads, 'payload_count': len(payloads),
      'payload_bytes': sum(v['bytes'] for v in payloads), 'deletion_recovery_ledger': deletions,
      'symlink_unlink_only': {'origin': str(resources.relative_to(ROOT)), 'target': os.readlink(resources)},
      'scope': 'Exact finished #94 scratch bytes and closure metadata recoverable from original trial/restoration or this local archive. Reviewed source also in pushed evidence tags. Ignored archives are local; no remote binary-backup or full qualification claim.'})
for v in deletions: assert seal(ROOT / v['origin']) == expected(v) == seal(ROOT / v['recovery'])
assert entries() == tree

def birth(pid):
    try: return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError): return None
known = []
owners = list(WORK.rglob('owner.json'))
assert len(owners) == 17
for path in owners:
    owner = read(path)
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    c = owner['controller']
    identities = [(c['pid'], c['start_ticks'])] + [(v['pid'], v['saved_start_ticks']) for v in owner['observed_births']]
    for pid, ticks in identities:
        current = birth(pid)
        assert current != ticks
        known.append({'pid': pid, 'saved_start_ticks': ticks, 'current_start_ticks': current, 'absent': True})
matches, errors, count = [], [], 0
for proc in Path('/proc').glob('[0-9]*'):
    if int(proc.name) == os.getpid(): continue
    count += 1
    for field in ('cmdline', 'environ', 'maps', 'cwd', 'exe'):
        try:
            value = os.readlink(proc / field) if field in ('cwd', 'exe') else (proc / field).read_bytes().decode(errors='replace')
            if str(WORK) in value or str(WORK.relative_to(ROOT)) in value:
                matches.append({'pid': int(proc.name), 'field': field})
        except OSError as e: errors.append({'pid': int(proc.name), 'field': field, 'error': e.errno})
    try:
        for fd in (proc / 'fd').iterdir():
            try:
                if str(WORK) in os.readlink(fd): matches.append({'pid': int(proc.name), 'field': 'fd', 'fd': fd.name})
            except OSError as e: errors.append({'pid': int(proc.name), 'field': 'fd/' + fd.name, 'error': e.errno})
    except OSError as e: errors.append({'pid': int(proc.name), 'field': 'fd', 'error': e.errno})
assert not matches, matches
check_retained()
assert entries() == tree
cleanup.mkdir(parents=True)
receipt = {'revision': head, 'source_committed_pushed': True, 'removed': False,
           'files': len(deletions), 'bytes': sum(v['bytes'] for v in deletions), 'targets': [str(WORK.relative_to(ROOT))],
           'publication_manifest': seal(publication / 'manifest.json'), 'CI_status_at_cleanup': ci['status'],
           'CI_conclusion_at_cleanup': ci['conclusion'], 'Linux_known_births': known,
           'Linux_accessible_processes': count, 'Linux_matches': matches, 'Linux_read_errors': errors,
           'resources_symlink_unlink_only': True, 'resource_target_files': 23,
           'scope': 'Only finished #94 scratch; resource symlink is unlinked without traversal. Ordinary builds/toolchains, archives, official exports and renders retained. Fresh finite known-birth/accessibility scan; denied fields, unobserved or scheduled/global consumers not proven. No Windows stage created by this batch.'}
write(cleanup / 'cleanup.json', receipt)
resources.unlink()
shutil.rmtree(WORK)
assert not WORK.exists() and target.is_dir()
check_retained()
receipt.update(removed=True, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               out_empty=not any((ROOT / 'out').iterdir()))
write(cleanup / 'cleanup.json', receipt)
print(json.dumps({'cleanup': str(cleanup.relative_to(ROOT)), 'files': receipt['files'], 'bytes': receipt['bytes'],
                  'out_empty': receipt['out_empty'], 'publication_manifest': receipt['publication_manifest']}))
