"""Protect finished rejected trial/restoration publications and remove only their staging."""
from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def read(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def identity(path):
    return {'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
def expected(record): return {k: record[k] for k in ('bytes', 'sha256')}
source = read(WORK / 'restoration/source.json')
head = source['source_revision']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert subprocess.check_output(['git', 'rev-parse', 'HEAD']) == subprocess.check_output(['git', 'rev-parse', '@{u}'])
for case in ('configure','build','admission'):
    owner=read(WORK/'restoration'/case/'owner.json')
    assert owner['passed'] and owner['child_exit']==0 and owner['source']=={'revision':head,'status':''}
assert read(WORK/'restoration/build-and-arrays.json')['original40_complete_arrays_and_whole_header_byte_exact95']
assert read(WORK/'restoration/linux-readonly-payloads.json')['pass_']
assert read(WORK/'native-preservation-reviewed.json')['accepted']
assert read(WORK/'software-preservation-reviewed.json')['accepted']
assert read(WORK/'restoration-preservation-reviewed.json')['accepted']
assert not read(WORK/'native-comparison.json')['pass']
ci=read(WORK/'restoration-ci.json')
assert ci['headSha']==head and ci['status'] in ('queued','in_progress','completed')
if ci['status']=='completed': assert ci['conclusion']=='success'
preserved=read(WORK/'tools-preservation-v3.json')
assert preserved['pushed'] and preserved['source_revision']==head
assert subprocess.check_output(['git','rev-parse',preserved['tag']+'^{commit}'],text=True).strip()==preserved['commit']
assert subprocess.check_output(['git','ls-remote','origin','refs/tags/'+preserved['tag']+'^{}'],text=True).split()[0]==preserved['commit']
for name,seal in preserved['files'].items():
    data=subprocess.check_output(['git','show',preserved['tag']+':'+name])
    assert {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}==expected(seal)
    assert identity(WORK/seal.get('origin',name))==expected(seal)
assert read(WORK / 'fresh-native-consumers-accepted.json')['pass']

recoveries = {}
native = ROOT / read(WORK / 'native-protection.json')['archive']
manifest = read(native / 'manifest.json')
assert identity(native / 'manifest.json') == read(WORK / 'native-protection.json')['manifest']
for record in (*manifest['executed_input_ledger'], *manifest['diagnostics']):
    recovery = record['recovery']
    if recovery['kind'] == 'immutable-Git-blob':
        data = subprocess.check_output(['git', 'cat-file', 'blob', recovery['blob']])
        assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected(record)
    else: assert identity(ROOT / recovery['path']) == expected(record)
    recoveries[record['origin']] = record
for name in ('software-protection.json','restoration-protection.json'):
    receipt = read(WORK / name)
    archive = ROOT / receipt['archive']
    assert identity(archive / 'manifest.json') == receipt['manifest']
    for record in read(archive / 'manifest.json')['payloads']:
        assert identity(archive / record['path']) == expected(record)
        recoveries[record['origin']] = {**record, 'recovery': {'kind': 'retained-archive-payload', 'path': str((archive / record['path']).relative_to(ROOT))}}

publication = ROOT / 'attestations/review' / head[:7] / 'rejected-fma-publication'
cleanup = ROOT / 'attestations/cleanup' / head[:7] / 'rejected-fma-finished'
assert not publication.exists() and not cleanup.exists()
publication_records = []
deletions = []
trial=read(WORK/'source.json')
assert trial['source_revision']==manifest['source_revision'] and trial['baseline_revision']==manifest['baseline_revision']
stages=[ROOT/'bin/windows-msvc'/('native-'+revision[:7]) for revision in (trial['baseline_revision'],trial['source_revision'])]
extras=[ROOT/'bin/linux-gcc/src/sirius/backend/retained'/name for name in ('retained_dense_fma.spv','retained_dopri_phase_fma.spv')]
directory_sets = {}
def entries(directory):
    values = {}
    assert directory.is_dir() and not directory.is_symlink()
    for path in directory.rglob('*'):
        assert not path.is_symlink()
        assert path.is_dir() or path.is_file()
        values[str(path.relative_to(directory))] = 'directory' if path.is_dir() else 'file'
    return values
for label, directory in [('observation', WORK), *[(p.name, p) for p in stages]]:
    assert directory.is_dir() and not directory.is_symlink()
    directory_sets[str(directory)] = entries(directory)
    for path in sorted(directory.rglob('*')):
        assert not path.is_symlink()
        if not path.is_file(): continue
        origin = str(path.relative_to(ROOT))
        value = identity(path)
        record = recoveries.get(origin)
        if record is not None and expected(record) == value:
            recovery = record['recovery']
        else:
            target = label + '/' + str(path.relative_to(directory))
            output = publication / target
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, output)
            assert identity(output) == value
            publication_records.append({'path': target, 'origin': origin, **value})
            recovery = {'kind': 'retained-publication-payload', 'path': str(output.relative_to(ROOT))}
        deletions.append({'origin': origin, **value, 'recovery': recovery})
for path in extras:
    assert path.is_file() and not path.is_symlink()
    origin=str(path.relative_to(ROOT));value=identity(path);record=recoveries[origin]
    assert expected(record)==value
    deletions.append({'origin':origin,**value,'recovery':record['recovery']})
publication.mkdir(parents=True, exist_ok=True)
(publication / 'manifest.json').write_text(json.dumps({'revision': head, 'payloads': publication_records,
    'payload_count': len(publication_records), 'payload_bytes': sum(v['bytes'] for v in publication_records),
    'deletion_recovery_ledger': deletions, 'scope': 'Exact final tool/review/CI/GitHub/observation metadata; references prior immutable raw and input archives. Local ignored evidence, not a remote backup or full qualification.'}, indent=2) + '\n')

targets = ['fma-dense-dp-review', *[p.name for p in stages]]
known = []
def birth(pid):
    try: return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError): return None
for directory in (WORK,):
    for path in directory.rglob('owner.json'):
        owner = read(path)
        if 'observed_births' not in owner: continue
        assert not owner['cleanup_errors'] and owner['remaining_owned_processes'] == []
        controller = owner['controller']
        current = birth(controller['pid'])
        assert current != controller['start_ticks']
        known.append({'pid': controller['pid'], 'saved_start_ticks': controller['start_ticks'], 'current_start_ticks': current, 'absent': True})
        for value in owner['observed_births']:
            current = birth(value['pid'])
            assert current != value['saved_start_ticks']
            known.append({'pid': value['pid'], 'saved_start_ticks': value['saved_start_ticks'], 'current_start_ticks': current, 'absent': True})
    for path in directory.glob('*-bridge.json'):
        owner = read(path)
        for pid_key, ticks_key in (('linux_owner_pid', 'linux_owner_start_ticks'), ('linux_bridge_child_pid', 'linux_bridge_child_start_ticks'), ('linux_child_pid', 'linux_child_birth')):
            if pid_key not in owner: continue
            current = birth(owner[pid_key])
            assert current != owner[ticks_key]
            known.append({'pid': owner[pid_key], 'saved_start_ticks': owner[ticks_key], 'current_start_ticks': current, 'absent': True})
for path in WORK.glob('native-*-sequence.json'):
    sequence=read(path)
    assert sequence['completed']
    checks=[(sequence['linux_owner_pid'],sequence['linux_owner_start_ticks'])]
    for call in sequence['owned_calls']:
        assert call['returncode']==0 and call['ownership_verified']
        checks.append((call['pid'],call['start_ticks']))
    for pid,ticks in checks:
        current=birth(pid);assert current!=ticks
        known.append({'pid':pid,'saved_start_ticks':ticks,'current_start_ticks':current,'absent':True})
windows = read(WORK / 'fresh-native-consumers.json')
assert windows['pass'] and not windows['matches'] and len(windows['known_births']) == 36
assert all(v['owned_birth_absent'] for v in windows['known_births'])
assert set(windows['targets']) == set(targets)
def check_fresh_windows():
    age = datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(windows['checked_utc'].replace('Z', '+00:00'))
    assert datetime.timedelta() <= age <= datetime.timedelta(seconds=60), 'Refresh Windows scoped census immediately before cleanup'
check_fresh_windows()
matches, errors, count = [], [], 0
for path in Path('/proc').glob('[0-9]*'):
    if int(path.name) == os.getpid(): continue
    count += 1
    for field in ('cmdline', 'environ', 'maps', 'cwd', 'exe'):
        try:
            value = os.readlink(path / field) if field in ('cwd', 'exe') else (path / field).read_bytes().decode(errors='replace')
            if any(target in value for target in [*targets,*[str(p.relative_to(ROOT)) for p in extras]]): matches.append({'pid': int(path.name), 'field': field})
        except (OSError, ValueError): errors.append({'pid': int(path.name), 'field': field})
assert not matches, matches
cleanup.mkdir(parents=True)
receipt = {'revision': head, 'source_committed_pushed': True, 'publication_manifest': identity(publication / 'manifest.json'),
           'removed': False, 'files': len(deletions), 'bytes': sum(v['bytes'] for v in deletions),
           'targets': [str(p.relative_to(ROOT)) for p in (WORK, *stages,*extras)],'CI_status_at_cleanup':ci['status'],'CI_conclusion_at_cleanup':ci['conclusion'],
           'Linux_known_births': known, 'Linux_accessible_processes': count, 'Linux_matches': matches, 'Linux_read_errors': errors,
           'Windows_scoped_census': windows,
           'scope': 'Finished rejected FMA trial/restoration scratch, its two obsolete native execution stages and two no-longer-generated shader byproducts only. Ordinary build/toolchain trees, official exports, raw archives and useful renders retained. Finite accessible scans; denied fields, handles and global/scheduled consumers unverified.'}
(cleanup / 'cleanup.json').write_text(json.dumps(receipt, indent=2) + '\n')
for record in deletions:
    assert identity(ROOT / record['origin']) == expected(record)
    recovery = record['recovery']
    if recovery['kind'] != 'immutable-Git-blob': assert identity(ROOT / recovery['path']) == expected(record)
for directory in (WORK, *stages): assert entries(directory) == directory_sets[str(directory)]
check_fresh_windows()
for path in extras: path.unlink()
for directory in (WORK, *stages): shutil.rmtree(directory)
receipt.update(removed=True, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), out_empty=not any((ROOT / 'out').iterdir()))
(cleanup / 'cleanup.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'cleanup': str(cleanup.relative_to(ROOT)), 'files': receipt['files'], 'bytes': receipt['bytes'], 'out_empty': receipt['out_empty'], 'publication_manifest': receipt['publication_manifest']}))
