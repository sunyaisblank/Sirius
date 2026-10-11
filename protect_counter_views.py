"""Push useful diagnostic source, then preserve complete finished capture inputs."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess, zipfile
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))
def seal(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
def expected(value):
    return {k: value[k] for k in ('bytes', 'sha256')}
def key(value):
    return value['bytes'], value['sha256']
def git(*args, **kwargs):
    return subprocess.check_output(['git', *args], **kwargs).decode().strip()
head = read(WORK / 'source.json')['source_revision']
assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}')
assert not git('status', '--porcelain')
assert git('rev-parse', 'HEAD^{tree}') == git('rev-parse', 'b81e061^{tree}')
assert len(sys.argv) == 2 and sys.argv[1] in ('sources', 'archive')
tag = 'evidence/issue-44-native-counter-views-' + head[:7]
DEST = ROOT / 'attestations/native-vulkan' / head[:7] / 'native-counter-views'
if sys.argv[1] == 'sources':
    review = read(WORK / 'preservation-source-review.json')
    assert review['accepted']
    for name, value in review['controller_seals'].items():
        assert seal(WORK / name) == value
    disposition = read(WORK / 'batch-disposition.json')
    assert disposition['finished_capture_batch'] and not disposition['full_qualification']
    assert not git('ls-remote', 'origin', 'refs/tags/' + tag)
    assert subprocess.run(['git', 'show-ref', '--verify', '--quiet',
                           'refs/tags/' + tag]).returncode != 0
    paths = sorted(p for p in WORK.iterdir() if p.is_file() and
                   (p.suffix in ('.py','.ps1','.cs','.cpp','.h','.md','.txt','.json') or
                    p.name in ('libamdrdf-license',)))
    index = WORK / 'preservation.index'
    assert not index.exists() and not (WORK / 'source-preservation.json').exists()
    env = dict(os.environ, GIT_INDEX_FILE=str(index))
    subprocess.run(['git', 'read-tree', '--empty'], env=env, check=True)
    records = []
    for path in paths:
        assert not path.is_symlink()
        data = path.read_bytes()
        blob = git('hash-object', '-w', '--stdin', input=data)
        name = str(path.relative_to(WORK))
        subprocess.run(['git', 'update-index', '--add', '--cacheinfo',
                        '100644,' + blob + ',' + name], env=env, check=True)
        records.append({'origin': str(path.relative_to(ROOT)), 'path': name,
                        'blob': blob, 'bytes': len(data),
                        'sha256': hashlib.sha256(data).hexdigest()})
    tree = git('write-tree', env=env)
    message = ('Preserve reviewed saved-profile counter and event observations at '
               + head[:7] + ' (#44)\n\n' + disposition['git_disposition'] + '\n')
    commit = git('commit-tree', tree, input=message.encode())
    subprocess.run(['git', 'tag', tag, commit], check=True)
    subprocess.run(['git', 'push', 'origin', 'refs/tags/' + tag], check=True)
    assert git('ls-remote', 'origin', 'refs/tags/' + tag).split()[0] == commit
    index.unlink()
    assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}')
    assert not git('status', '--porcelain')
    receipt = {'tag': tag, 'commit': commit, 'source_revision': head,
               'pushed': True, 'remote_commit_verified': True, 'files': records}
    (WORK / 'source-preservation.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'tag': tag, 'commit': commit, 'source_files': len(records),
                      'pushed': True}))
    sys.exit(0)

source = read(WORK / 'source-preservation.json')
assert source['pushed'] and source['remote_commit_verified'] and source['tag'] == tag
assert git('rev-parse', tag) == source['commit']
assert git('ls-remote', 'origin', 'refs/tags/' + tag).split()[0] == source['commit']
assert not DEST.exists()
recoveries = {}
prior_path = ROOT / 'attestations/native-vulkan/6be3132/native-spm-cost/manifest.json'
prior = read(prior_path)
for value in [*prior['executed_input_ledger'], *prior['records'], prior['tool_archive']]:
    recoveries[key(value)] = value['recovery']
recoveries[key(seal(prior_path))] = {'kind':'retained-prior-payload','path':str(prior_path.relative_to(ROOT))}
for value in source['files']:
    recoveries[key(value)] = {'kind': 'immutable-Git-blob',
                             'revision': source['commit'], 'blob': value['blob'],
                             'path': value['path']}
payloads = {}
verified = set()
def verify_recovery(recovery, value):
    identity = json.dumps(recovery, sort_keys=True), key(value)
    if identity in verified:
        return
    if recovery['kind'] == 'immutable-Git-blob':
        data = subprocess.check_output(['git', 'cat-file', 'blob', recovery['blob']])
        assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected(value)
    elif recovery['kind'] == 'protected-ZIP-member':
        archive = ROOT / recovery['archive_path']
        assert seal(archive) == recovery['archive_seal']
        with zipfile.ZipFile(archive) as z:
            info = z.getinfo(recovery['member'])
            assert info.file_size == value['bytes'] and info.CRC == recovery['CRC32']
            data = z.read(info)
            assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected(value)
    else:
        assert seal(ROOT / recovery['path']) == expected(value)
    verified.add(identity)
def resolve(path, value):
    identity = key(value)
    if identity not in recoveries:
        assert path.is_file() and not path.is_symlink() and seal(path) == expected(value), str(path)
        relative = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else None
        tracked = git('ls-files', '--', relative) if relative else ''
        if tracked == relative:
            blob = git('rev-parse', head + ':' + relative)
            recoveries[identity] = {'kind': 'immutable-Git-blob', 'revision': head,
                                    'blob': blob, 'path': relative}
        else:
            target = DEST / 'payloads' / value['sha256']
            target.parent.mkdir(parents=True, exist_ok=True)
            assert not target.exists()
            shutil.copy2(path, target)
            assert seal(target) == expected(value)
            payloads[str(target.relative_to(DEST))] = expected(value)
            recoveries[identity] = {'kind': 'local-archive-payload',
                                    'path': str(target.relative_to(ROOT))}
    recovery = recoveries[identity]
    verify_recovery(recovery, value)
    return recovery

download = prior['tool_archive']
assert download['recovery']['kind'] == 'local-archive-payload'
archive = ROOT / download['recovery']['path']
archive_value = expected(download)
assert seal(archive) == archive_value
archive_recovery = resolve(archive, archive_value)
protected_archive = archive_recovery['path']
tools = [*read(WORK / 'tools-extraction.json')['files'],
         *read(WORK / 'tools-gui-extraction.json')['files']]
assert len(tools) == 70 and len({v['path'] for v in tools}) == 70
with zipfile.ZipFile(ROOT / protected_archive) as z:
    for value in tools:
        assert seal(ROOT / value['path']) == expected(value)
        info = z.getinfo(value['zip_member'])
        assert info.file_size == value['bytes'] and info.CRC == value['CRC32']
        data = z.read(info)
        assert {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} == expected(value)
        recovery = {'kind': 'protected-ZIP-member', 'archive_path': protected_archive,
                    'archive_seal': archive_value, 'member': value['zip_member'],
                    'CRC32': value['CRC32']}
        recoveries[key(value)] = recovery
        verified.add((json.dumps(recovery, sort_keys=True), key(value)))

ledger = {}
def collect(origin, value, provenance):
    identity = origin, *key(value)
    path = Path(origin) if Path(origin).is_absolute() else ROOT / origin
    current = path.is_file() and seal(path) == expected(value)
    recovery = resolve(path, value)
    if identity not in ledger:
        ledger[identity] = {'origin': origin, **expected(value),
                            'origin_current_at_preservation': current,
                            'ledgers': [], 'recovery': recovery}
    ledger[identity]['ledgers'].append(provenance)
for path in sorted(WORK.rglob('*.json')):
    value = read(path)
    def visit(node):
        if isinstance(node, dict):
            for name, child in node.items():
                if name in ('input_seals', 'source_input_seals', 'files') and isinstance(child, dict):
                    for origin, item in child.items():
                        if isinstance(item, dict) and set(('bytes', 'sha256')) <= item.keys():
                            collect(origin, item, str(path.relative_to(ROOT)))
                else:
                    visit(child)
        elif isinstance(node, list):
            for child in node:
                visit(child)
    visit(value)
records = []
targets = [WORK, ROOT / 'bin/windows-msvc/native-b81e061-executed-cost',
           ROOT / 'bin/windows-msvc/qualification-tools/rdts-2026-10-02-1872']
for target in targets:
    assert target.is_dir() and not target.is_symlink()
    for path in sorted(target.rglob('*')):
        assert not path.is_symlink(), str(path)
        if path.is_file():
            value = seal(path)
            records.append({'origin': str(path.relative_to(ROOT)), **value,
                            'recovery': resolve(path, value)})
stage = [v for v in records if v['origin'].startswith('bin/windows-msvc/native-b81e061-executed-cost/')]
assert len(stage) == 94
assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}') and not git('status', '--porcelain')
manifest = {'source_revision': head, 'baseline_revision': 'b81e061fa081cdaed37eca3add950bd254bf296d',
            'source_snapshot': source, 'disposition': read(WORK / 'batch-disposition.json'),
            'prior_manifest': {'path': str(prior_path.relative_to(ROOT)), **seal(prior_path)},
            'executed_input_ledger': list(ledger.values()), 'records': records,
            'tool_archive': {'origin': str(archive.relative_to(ROOT)), **archive_value,
                             'recovery': archive_recovery},
            'tools': tools, 'payloads': payloads,
            'payload_count': len(payloads), 'payload_bytes': sum(v['bytes'] for v in payloads.values()),
            'scope': 'Whole finished saved-profile counter diagnostic files and exact consumed input versions preserved. Prior acquisition trace remains immutable. Useful source is pushed; raw observation evidence and existing official tool ZIP remain local. No ordinary performance benefit, rendering/scientific/release completion or cleanup clearance.'}
DEST.mkdir(parents=True, exist_ok=True)
(DEST / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
receipt = {'archive': str(DEST.relative_to(ROOT)), 'manifest': seal(DEST / 'manifest.json'),
           'records': len(records), 'executed_origin_version_pairs': len(ledger),
           'payload_count': len(payloads), 'payload_bytes': manifest['payload_bytes'],
           'all_recovery_bytes_verified': True, 'cleanup_authorized': False}
(WORK / 'protection.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt))
