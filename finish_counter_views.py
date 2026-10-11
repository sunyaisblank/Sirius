"""Verify recovery bytes, preserve final receipts, and remove only finished scratch."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess, zipfile
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
assert ROOT == Path('/home/astra/.project/Sirius')
assert WORK == ROOT / 'out/native-counter-views'
assert len(sys.argv) == 2 and sys.argv[1] in ('review', 'remove')
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def expected(v): return {k: v[k] for k in ('bytes', 'sha256')}
def key(v): return v['bytes'], v['sha256']
def git(*a): return subprocess.check_output(['git', *a], text=True).strip()
def birth(pid):
    try: return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError): return None
head = read(WORK / 'source.json')['source_revision']
def source_check():
    assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}')
    assert not git('status', '--porcelain')
    assert git('rev-parse', 'HEAD^{tree}') == git('rev-parse', 'b81e061^{tree}')
    assert git('rev-parse', 'main') == git('rev-parse', 'origin/main') == '5e706073a3fd2f9afcf8efbcde2f094122050f6a'
    assert len([v for v in git('worktree', 'list', '--porcelain').splitlines() if v.startswith('worktree ')]) == 1
source_check()
review = read(WORK / 'preservation-source-review.json'); assert review['accepted']
for name, value in review['controller_seals'].items(): assert seal(WORK / name) == value
protection = read(WORK / 'protection.json')
DEST = ROOT / protection['archive']
assert DEST == ROOT / 'attestations/native-vulkan' / head[:7] / 'native-counter-views'
assert seal(DEST / 'manifest.json') == protection['manifest']
manifest = read(DEST / 'manifest.json')
source = read(WORK / 'source-preservation.json')
assert source == manifest['source_snapshot'] and source['pushed'] and source['remote_commit_verified']
assert git('rev-parse', source['tag']) == source['commit']
assert git('ls-remote', 'origin', 'refs/tags/' + source['tag']).split()[0] == source['commit']
verified, archives, recoveries = set(), {}, {}
def verify(v, r):
    marker = key(v), json.dumps(r, sort_keys=True)
    if marker in verified: return
    kind = r['kind']
    if kind == 'immutable-Git-blob':
        assert git('rev-parse', r['revision'] + ':' + r['path']) == r['blob']
        data = subprocess.check_output(['git', 'cat-file', 'blob', r['blob']])
        assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected(v)
    elif kind == 'protected-ZIP-member':
        p = ROOT / r['archive_path']
        if str(p) not in archives:
            assert p.is_file() and not p.is_symlink() and seal(p) == r['archive_seal']
            archives[str(p)] = zipfile.ZipFile(p)
        z = archives[str(p)]; info = z.getinfo(r['member'])
        assert info.file_size == v['bytes'] and info.CRC == r['CRC32']
        data = z.read(info)
        assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected(v)
    else:
        assert kind in ('local-archive-payload', 'local-publication-payload', 'retained-native-export', 'retained-prior-payload', 'retained-provider-payload', 'retained-archive-payload'), kind
        p = ROOT / r['path']; assert p.is_file() and not p.is_symlink() and seal(p) == expected(v)
    verified.add(marker); recoveries[key(v)] = r
for v in [*manifest['executed_input_ledger'], *manifest['records'], manifest['tool_archive']]: verify(v, v['recovery'])
for name, value in manifest['payloads'].items(): assert seal(DEST / name) == value
for v in source['files']:
    assert seal(ROOT / v['origin']) == expected(v)
    verify(v, {'kind':'immutable-Git-blob', 'revision':source['commit'], 'blob':v['blob'], 'path':v['path']})
assert len([v for v in manifest['records'] if v['origin'].startswith('bin/windows-msvc/native-b81e061-executed-cost/')]) == 94
assert len(manifest['tools']) == 70
for z in archives.values(): z.close()
if sys.argv[1] == 'review':
    p = WORK / 'preservation-reviewed.json'; assert not p.exists()
    p.write_text(json.dumps({'accepted':True, 'manifest':protection['manifest'], 'recovery_versions':len(verified), 'source_tag':source['tag'], 'pushed_source_verified':True, 'scope':'Independent full-byte Git/payload/ZIP-member recovery verification; no cleanup or qualification verdict.'}, indent=2)+'\n')
    print(json.dumps({'accepted':True, 'recovery_versions':len(verified)})); sys.exit(0)
assert read(WORK / 'preservation-reviewed.json')['accepted']
assert read(WORK / 'preservation-reviewed.json')['manifest'] == protection['manifest']
census = read(WORK / 'fresh-finished-consumers-bridge.json')
windows = read(WORK / 'fresh-finished-consumers/owner.json')
assert census['accepted'] and census['returncode'] == 0 and not census['errors'] and not census['outer_stop'] and not census['controller_signals']
assert windows['pass'] and not windows['errors'] and not windows['matches']
known = read(WORK / 'finished-births.json')
assert sorted((v['pid'],v['birth']) for v in windows['known_births']) == sorted((v['pid'],v['birth']) for v in known['windows'])
assert all(v['owned_birth_absent'] for v in windows['known_births'])
emergency = read(WORK / 'fresh-finished-consumers-emergency.json')
assert emergency['recorded_births_absent'] and not emergency['errors'] and not any(v['same_birth_live'] for v in emergency['targets'])
def fresh():
    age = datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(windows['checked_utc'].replace('Z','+00:00'))
    assert datetime.timedelta() <= age <= datetime.timedelta(seconds=60), age
fresh()
linux = []
for v in [*known['linux'], {'pid':census['linux_owner_pid'],'start_ticks':census['linux_owner_start_ticks']}, {'pid':census['linux_child_pid'],'start_ticks':census['linux_child_start_ticks']}, census['emergency_observer']]:
    ticks = str(v['start_ticks']); assert birth(v['pid']) != ticks, v
    linux.append({'pid':v['pid'],'start_ticks':ticks,'absent':True})
bootstrap = read(WORK / 'fresh-finished-consumers-cleanup-helper-bootstrap.json')
assert bootstrap['owner_pid'] == emergency['helper_pid'] and bootstrap['owner_start_time_utc'] == emergency['helper_start_time_utc']
# Emergency observer returned normally; no extra leaf observer was needed.
assert 'cleanup_stop_observer' not in census
plan = manifest['disposition']['cleanup']
zip_origin = manifest['tool_archive']['origin']
assert plan['remove'] == ['out/native-counter-views']
assert zip_origin == 'attestations/native-vulkan/6be3132/native-executed-cost/payloads/128fdebb889d1fda1d04ba9e4792fbef1b18c086b4e901642262653a7a53917e'
assert {v['path'] for v in plan['retain']} == {'bin/windows-msvc/native-b81e061-executed-cost','bin/windows-msvc/qualification-tools/rdts-2026-10-02-1872'}
def tree(folder):
    assert folder.is_dir() and not folder.is_symlink()
    result = {}
    for p in sorted(folder.rglob('*')):
        assert not p.is_symlink() and (p.is_dir() or p.is_file())
        result[str(p.relative_to(ROOT))] = seal(p) if p.is_file() else 'directory'
    return result
retained = {v['origin']: expected(v) for v in manifest['records'] if v['origin'].startswith('bin/') and v['origin'] != zip_origin}
assert len(retained) == 164
for name, value in retained.items(): assert seal(ROOT / name) == value
resources = {str(p): tree(ROOT / p) for p in (Path('bin/linux-gcc/src/sirius/backend/retained'),Path('bin/linux-gcc/tests/backend/resources'),Path('bin/linux-gcc/src/sirius/app/resources'))}
assert [sum(isinstance(v,dict) for v in values.values()) for values in resources.values()] == [35,23,32]
closing = DEST / 'closing'; assert not closing.exists()
original_tree = tree(WORK)
entries = []
for p in [ROOT / name for name,v in original_tree.items() if isinstance(v,dict)]:
    assert p.is_file() and not p.is_symlink(); value = seal(p); r = recoveries.get(key(value))
    if r is None:
        q = closing / 'payloads' / value['sha256']; q.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,q)
        assert seal(q) == value
        r = {'kind':'local-publication-payload','path':str(q.relative_to(ROOT))}
    verify(value,r)
    entries.append({'origin':str(p.relative_to(ROOT)),**value,'recovery':r})
assert tree(WORK) == original_tree
matches, errors, count = [], [], 0
for p in Path('/proc').glob('[0-9]*'):
    if int(p.name) == os.getpid(): continue
    count += 1
    for field in ('cmdline','environ','maps','cwd','exe'):
        try:
            value = os.readlink(p / field) if field in ('cwd','exe') else (p/field).read_bytes().decode(errors='replace')
            if any(v in value for v in ('native-counter-views', 'native-b81e061-executed-cost', 'rdts-2026-10-02-1872')): matches.append({'pid':int(p.name),'field':field})
        except (OSError,ValueError): errors.append({'pid':int(p.name),'field':field})
assert not matches, matches
source_check(); fresh()
closing.mkdir(parents=True,exist_ok=True)
closing_manifest = {'source_revision':head, 'source_tag':source, 'deletion_recovery_ledger':entries, 'retained_products_and_tools':retained, 'canonical_product_resource_trees':resources, 'scope':'Finished profiler diagnostic closure; pushed useful source and exact local raw recovery. No new optimization, frame or scientific/release qualification.'}
(closing / 'manifest.json').write_text(json.dumps(closing_manifest,indent=2)+'\n')
receipt = {'source_revision':head, 'source_committed_pushed':True, 'closing_manifest':seal(closing/'manifest.json'), 'files':len(entries), 'bytes':sum(v['bytes'] for v in entries), 'remove':plan['remove'], 'retain':plan['retain'], 'linux_recorded_births':linux, 'linux_accessible_processes':count,'linux_matches':matches,'linux_read_errors':errors,'windows_census':windows,'removed':False,'scope':'Only finished saved-profile counter diagnostic scratch removed; existing protected ZIP preserved. Retained 94 native files, 70 extracted official tools and canonical products/resources. Finite accessible fields only; hidden/global/scheduled consumers remain unverified.'}
(DEST / 'cleanup.json').write_text(json.dumps(receipt,indent=2)+'\n')
assert tree(WORK) == original_tree
assert seal(ROOT / zip_origin) == expected(manifest['tool_archive'])
for name,value in retained.items(): assert seal(ROOT/name) == value
for name,value in resources.items(): assert tree(ROOT/name) == value
source_check(); fresh()
shutil.rmtree(WORK)
assert not WORK.exists() and seal(ROOT / zip_origin) == expected(manifest['tool_archive'])
for name,value in retained.items(): assert seal(ROOT/name) == value
for name,value in resources.items(): assert tree(ROOT/name) == value
source_check()
receipt.update(removed=True, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),retained_verified_after=True,canonical_verified_after=True,out_empty=not any((ROOT/'out').iterdir()),final_tracking=git('rev-list','--left-right','--count','HEAD...@{u}'))
(DEST/'cleanup.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('files','bytes','removed','out_empty','final_tracking')}))
