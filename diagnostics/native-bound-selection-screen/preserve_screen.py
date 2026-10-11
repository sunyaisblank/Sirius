"""Preserve useful source and raw screen outputs before finished WORK removal."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HEAD = '2053ef79c90752d37b0075b8137cfb9b571b6820'
ARCHIVE = ROOT / 'attestations/native-vulkan/2053ef7/native-bound-selection-screen'
TAG = 'evidence/issue-44-native-bound-selection-screen-2053ef7'
PREFIX = 'diagnostics/native-bound-selection-screen/'

def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

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
    assert WORK == ROOT / 'out/native-bound-selection' and not WORK.is_symlink()
    assert git('rev-parse', 'HEAD') == HEAD == git('rev-parse', '@{u}')
    assert git('branch', '--show-current') == 'development/retained-dense-step-control-2026-10-06'
    assert not git('status', '--porcelain')
    assert git('rev-list', '--left-right', '--count', 'main...origin/main') == '0\t0'
    assert git('worktree', 'list', '--porcelain').count('worktree ') == 1
    assert read(WORK / 'findings.json')['accepted']
    assert not read(WORK / 'findings.json')['source_optimization_adopted']
    for directory in (WORK, WORK / 'select'):
        for name in ('source-review.json', 'compile-admission.json', 'compiled-review.json',
                     'driver-review.json', 'ownership-review.json'):
            assert read(directory / name)['accepted']
        compiler = read(directory / 'compiler.json')
        assert len(compiler['runs']) == 4 and all(run['accepted'] for run in compiler['runs'])
        assert not compiler['signals']
        for name, expected in compiler['input_seals'].items():
            assert seal(Path(name)) == expected
        admission = read(directory / 'compile-admission.json')
        for name, key in (('phase/retained_transport.spv', 'spv'),
                          ('phase/retained_transport.spvasm', 'assembly'),
                          ('compiler.json', 'compiler_receipt')):
            assert seal(directory / name) == admission[key]
        assert seal(directory / 'compile-admission.json') == read(directory / 'compiled-review.json')['admission_seal']
        assert seal(directory / 'inspect_compiled.py') == read(directory / 'compiled-review.json')['inspector_seal']
        driver_review = read(directory / 'driver-review.json')
        assert seal(directory / 'driver-inspection.json') == driver_review['inspection_seal']
        assert seal(directory / 'inspect_driver.py') == driver_review['inspector_seal']
        for name, expected in read(directory / 'ownership-review.json')['input_seals'].items():
            assert seal(directory / name) == expected
        placement = directory / 'placement'
        bridge = read(placement / 'placement-query-bridge.json')
        assert bridge['accepted'] and bridge['returncode'] == 0
        for name, expected in bridge['input_seals'].items():
            assert seal(Path(name) if Path(name).is_absolute() else ROOT / name) == expected
        for name, expected in bridge['outputs'].items():
            assert seal(Path(name) if Path(name).is_absolute() else ROOT / name) == expected
        assert read(placement / 'retirement-run.json')['accepted']
        assert read(placement / 'native-retirement.json')['pass']

assert len(sys.argv) == 2 and sys.argv[1] in ('archive', 'remove')
guard()
if sys.argv[1] == 'archive':
    review = read(WORK / 'preservation-source-review.json')
    assert review['accepted'] and seal(Path(__file__)) == review['helper_seal']
    assert not ARCHIVE.exists() and not git('tag', '--list', TAG)
    env = dict(os.environ, GIT_INDEX_FILE=str(WORK / 'preservation.index'))
    blobs = {}
    useful_json = {'source-review.json', 'compiler-controller-review.json', 'compiled-review.json',
                   'driver-review.json', 'ownership-review.json', 'findings.json',
                   'preservation-source-review.json'}
    try:
        git('read-tree', HEAD, env=env)
        for path in files():
            relative = str(path.relative_to(WORK))
            useful = (path.suffix in ('.py', '.ps1', '.cs', '.patch') or
                      path.name in useful_json or relative in (
                          'kernels/retained_pair.slang', 'phase/retained_transport.slang',
                          'select/kernels/retained_pair.slang'))
            if 'native-temp' in path.parts or 'owner-temp' in path.parts:
                useful = False
            if not useful:
                continue
            blob = git('hash-object', '-w', str(path))
            name = PREFIX + relative
            git('update-index', '--add', '--cacheinfo', '100644', blob, name, env=env)
            blobs[relative] = {'git_path': name, 'git_blob': blob}
        tree = git('write-tree', env=env)
        commit = git('commit-tree', tree, '-p', HEAD, input=
                     'Preserve native bound selection compiler and driver screens (#44)\n\n'
                     'Retain both exact-word source candidates and actual static code evidence. '
                     'Production source remains unchanged. Two ordinary pipelines, zero dispatches; '
                     'numerical identity and complete-path benefit remain pending.\n')
        git('tag', TAG, commit)
    finally:
        for path in WORK.glob('preservation.index*'):
            path.unlink()
    ARCHIVE.mkdir(parents=True)
    records = []
    for path in files():
        target = ARCHIVE / 'files' / path.relative_to(WORK)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        assert seal(path) == seal(target)
        records.append({'origin': str(path.relative_to(ROOT)),
                        'recovery': str(target.relative_to(ROOT)), **seal(path),
                        **blobs.get(str(path.relative_to(WORK)), {})})
    manifest = {'source_revision': HEAD, 'source_tag': TAG, 'source_commit': commit,
                'records': records, 'scope': 'Complete production parent source plus useful '
                'diagnostic source/concise findings in Git; all raw outputs copied locally. '
                'No remote backup of raw native evidence or adopted optimization.'}
    (ARCHIVE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'source_tag': TAG, 'source_commit': commit,
                      'files': len(records), 'bytes': sum(row['bytes'] for row in records)}))
else:
    assert read(WORK / 'preservation-review.json')['accepted']
    assert read(WORK / 'publication-readback.json')['exact']
    final = read(WORK / 'final-retirement.json')
    assert final['accepted'] and final['linux_stage'] == str(WORK)
    assert seal(WORK / 'retirement-run.json') == final['bridge_seal']
    assert seal(WORK / 'native-retirement.json') == final['census_seal']
    final_bridge = read(WORK / 'retirement-run.json')
    assert final_bridge['accepted'] and not final_bridge['errors']
    assert final_bridge['arguments'][final_bridge['arguments'].index('-Stage') + 1] == final['windows_stage']
    expected_inputs = {str(WORK / name) for name in (
        'run_native_retirement.py', 'verify_native_retirement.ps1',
        'placement/run/owner.json', 'placement/owner-temp/bootstrap.json',
        'select/placement/run/owner.json', 'select/placement/owner-temp/bootstrap.json')}
    assert set(final_bridge['input_seals']) == expected_inputs
    for name, expected in final_bridge['input_seals'].items():
        assert seal(Path(name)) == expected
    census = read(WORK / 'native-retirement.json')
    assert census['pass'] and not census['errors'] and not census['matches']
    assert len(census['targets']) == 4 and all(target['recorded_birth_absent'] for target in census['targets'])
    manifest = read(ARCHIVE / 'manifest.json')
    assert git('ls-remote', 'origin', 'refs/tags/' + TAG).split()[0] == manifest['source_commit']
    for row in manifest['records']:
        assert seal(ROOT / row['origin']) == seal(ROOT / row['recovery']) == {
            key: row[key] for key in ('bytes', 'sha256')}
    records = []
    for path in files():
        target = ARCHIVE / 'files' / path.relative_to(WORK)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            assert seal(path) == seal(target)
        else:
            shutil.copyfile(path, target)
        assert seal(path) == seal(target)
        records.append({'origin': str(path.relative_to(ROOT)),
                        'recovery': str(target.relative_to(ROOT)), **seal(path)})
    for directory in Path('/proc').iterdir():
        if not directory.name.isdigit() or int(directory.name) == os.getpid():
            continue
        try:
            command = (directory / 'cmdline').read_bytes().replace(b'\0', b' ')
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        assert str(WORK).encode() not in command and b'out/native-bound-selection/' not in command
    publication = ARCHIVE / 'publication'
    publication.mkdir()
    ledger = {'source_tag': TAG, 'source_commit': manifest['source_commit'],
              'manifest': seal(ARCHIVE / 'manifest.json'), 'records': records}
    (publication / 'closing-ledger.json').write_text(json.dumps(ledger, indent=2) + '\n')
    shutil.rmtree(WORK)
    state = {'source_revision': HEAD, 'clean': not git('status', '--porcelain'),
             'tracking': git('rev-list', '--left-right', '--count', 'HEAD...@{u}'),
             'main_tracking': git('rev-list', '--left-right', '--count', 'main...origin/main'),
             'worktrees': git('worktree', 'list'), 'finished_work_absent': not WORK.exists(),
             'out_entries': sorted(p.name for p in (ROOT / 'out').iterdir()),
             'removed_files': len(records), 'removed_bytes': sum(row['bytes'] for row in records),
             'scope': 'Source/static driver screen complete; runtime and full qualification open.'}
    assert state['clean'] and state['tracking'] == state['main_tracking'] == '0\t0'
    (publication / 'removal-state.json').write_text(json.dumps(state, indent=2) + '\n')
    print(json.dumps(state))
