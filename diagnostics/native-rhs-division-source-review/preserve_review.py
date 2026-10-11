"""Preserve finite reviewed source findings, then remove only their finished scratch."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
ARCHIVE = ROOT / 'attestations/source-review/2053ef7/native-rhs-division-source-review'
INPUT_ARCHIVE = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution'
HEAD = '2053ef79c90752d37b0075b8137cfb9b571b6820'
TAG = 'evidence/issue-44-native-division-source-review-2053ef7'
PREFIX = 'diagnostics/native-rhs-division-source-review/'


def read(path):
    return json.loads(path.read_text())


def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def git(*args, env=None, input=None):
    return subprocess.check_output(['git', *args], env=env, input=input,
                                   text=True, timeout=30).strip()


def files():
    entries = list(WORK.rglob('*'))
    assert all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in entries)
    return sorted(p for p in entries if p.is_file())


def guard():
    assert WORK == ROOT / 'out/native-rhs-division-source-review' and not WORK.is_symlink()
    assert git('rev-parse', 'HEAD') == HEAD == git('rev-parse', '@{u}')
    assert git('branch', '--show-current') == 'development/retained-dense-step-control-2026-10-06'
    assert not git('status', '--porcelain')
    assert git('rev-list', '--left-right', '--count', 'main...origin/main') == '0\t0'
    assert git('worktree', 'list', '--porcelain').count('worktree ') == 1
    for name in ('findings.json', 'metadata-review.json', 'historical-review.json'):
        assert read(WORK / name)['accepted']
    metadata = read(WORK / 'metadata-review.json')
    assert seal(WORK / 'inspect_native_division.py')['sha256'] == metadata['reader_sha256']
    assert seal(WORK / 'native-division-inspection.json')['sha256'] == metadata['inspection_sha256']
    assert seal(INPUT_ARCHIVE / 'manifest.json')['sha256'] == \
        '71e2994b331a5a1e8aef3e3e5e472e19b1ce72e9b9c2aaed7430efff08510d2d'
    inputs = read(WORK / 'native-division-inspection.json')['inputs']
    assert set(inputs) == {'phase/retained_transport.spvasm',
                           'placement/query-payloads/transport.info1',
                           'placement/query-payloads/transport.info2'}
    for name, expected in inputs.items():
        assert seal(INPUT_ARCHIVE / 'files' / name) == expected
    for expected in read(WORK / 'historical-review.json')['evidence'].values():
        assert seal(ROOT / expected['path']) == {k: expected[k] for k in ('bytes', 'sha256')}
    assert not read(WORK / 'findings.json')['source_optimization_adopted']


assert len(sys.argv) == 2 and sys.argv[1] in ('archive', 'remove')
guard()
if sys.argv[1] == 'archive':
    assert not ARCHIVE.exists() and not git('tag', '--list', TAG)
    env = dict(os.environ, GIT_INDEX_FILE=str(WORK / 'preservation.index'))
    blobs = {}
    try:
        git('read-tree', HEAD, env=env)
        useful_names = {'findings.json', 'metadata-review.json', 'historical-review.json',
                        'preservation-source-review.json'}
        for path in files():
            if path.suffix != '.py' and path.name not in useful_names:
                continue
            name = PREFIX + str(path.relative_to(WORK))
            blob = git('hash-object', '-w', str(path))
            git('update-index', '--add', '--cacheinfo', '100644', blob, name, env=env)
            blobs[str(path.relative_to(WORK))] = {'git_path': name, 'git_blob': blob}
        tree = git('write-tree', env=env)
        commit = git('commit-tree', tree, '-p', HEAD, input=
                     'Preserve native division source attribution and closed prior leads (#44)\n\n'
                     'Join current diagnostic IR/ISA to exact ELF bytes and branch offsets. '
                     'Denominator absolute/exponent/CLZ are retained across the full correction loop; '
                     'eleven invariant integer steps remain without a material benefit premise. '
                     'NativeTransport32 and zero-quotient source screens are already closed. '
                     'No candidate, native run or scientific/release qualification.\n')
        git('tag', TAG, commit)
    finally:
        for path in WORK.glob('preservation.index*'):
            path.unlink()
    ARCHIVE.mkdir(parents=True)
    rows = []
    for path in files():
        target = ARCHIVE / 'files' / path.relative_to(WORK)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        assert seal(path) == seal(target)
        rows.append({'origin': str(path.relative_to(ROOT)),
                     'recovery': str(target.relative_to(ROOT)), **seal(path),
                     **blobs.get(str(path.relative_to(WORK)), {})})
    manifest = {'source_revision': HEAD, 'source_tag': TAG, 'source_commit': commit,
                'records': rows, 'scope': 'Complete finite source-review WORK recovery; '
                'complete production source plus useful reader and concise findings in Git. '
                'Raw compiled-code inspection and prior native evidence remain local.'}
    (ARCHIVE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'source_tag': TAG, 'source_commit': commit,
                      'files': len(rows), 'bytes': sum(r['bytes'] for r in rows)}))
else:
    assert read(WORK / 'preservation-review.json')['accepted']
    assert read(WORK / 'publication-readback.json')['exact']
    manifest = read(ARCHIVE / 'manifest.json')
    assert git('ls-remote', 'origin', 'refs/tags/' + TAG).split()[0] == manifest['source_commit']
    for row in manifest['records']:
        assert seal(ROOT / row['origin']) == seal(ROOT / row['recovery']) == \
            {k: row[k] for k in ('bytes', 'sha256')}
    rows = []
    for path in files():
        target = ARCHIVE / 'files' / path.relative_to(WORK)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            assert seal(path) == seal(target)
        else:
            shutil.copyfile(path, target)
        assert seal(path) == seal(target)
        rows.append({'origin': str(path.relative_to(ROOT)),
                     'recovery': str(target.relative_to(ROOT)), **seal(path)})
    # This batch launched no native/build processes. Check accessible Linux
    # consumers before WORK-only removal, without a global/scheduled claim.
    for directory in Path('/proc').iterdir():
        if not directory.name.isdigit() or int(directory.name) == os.getpid():
            continue
        try:
            cmd = (directory / 'cmdline').read_bytes().replace(b'\0', b' ')
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        assert str(WORK).encode() not in cmd and b'out/native-rhs-division-source-review/' not in cmd
    publication = ARCHIVE / 'publication'
    publication.mkdir()
    ledger = {'source_tag': TAG, 'source_commit': manifest['source_commit'],
              'manifest': seal(ARCHIVE / 'manifest.json'), 'records': rows}
    (publication / 'closing-ledger.json').write_text(json.dumps(ledger, indent=2) + '\n')
    shutil.rmtree(WORK)
    state = {'clean': not git('status', '--porcelain'), 'source_revision': HEAD,
             'tracking': git('rev-list', '--left-right', '--count', 'HEAD...@{u}'),
             'main_tracking': git('rev-list', '--left-right', '--count', 'main...origin/main'),
             'worktrees': git('worktree', 'list'), 'finished_work_absent': not WORK.exists(),
             'out_entries': sorted(p.name for p in (ROOT / 'out').iterdir()),
             'removed_files': len(rows), 'removed_bytes': sum(r['bytes'] for r in rows),
             'scope': 'Source review closed; long-term rendering and qualification remain open.'}
    assert state['clean'] and state['tracking'] == state['main_tracking'] == '0\t0'
    (publication / 'removal-state.json').write_text(json.dumps(state, indent=2) + '\n')
    print(json.dumps(state))
