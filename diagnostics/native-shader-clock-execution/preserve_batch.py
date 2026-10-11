"""Preserve reviewed diagnostic source and all finite raw evidence before WORK removal."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
ARCHIVE = ROOT / 'attestations/native-vulkan/8f1118d/native-shader-clock-execution'
HEAD = '8f1118dd56d10cdd0f932613427b68409a1f77f6'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
TAG = 'evidence/issue-44-native-shader-clock-execution-8f1118d'
PREFIX = 'diagnostics/native-shader-clock-execution/'
CASES = ('device-timestamps', 'mixed-transport', 'factored-endpoints')
def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*args, env=None, input=None):
    return subprocess.check_output(['git', *args], env=env, input=input, text=True, timeout=30).strip()
def birth(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None
def files():
    entries = list(WORK.rglob('*'))
    assert all(not p.is_symlink() for p in entries)
    return sorted(p for p in entries if p.is_file())
def guard():
    assert WORK == ROOT / 'out/native-shader-clock-execution' and not WORK.is_symlink()
    assert git('rev-parse', 'HEAD') == HEAD == git('rev-parse', '@{u}')
    assert not git('status', '--porcelain')
    assert git('rev-list', '--left-right', '--count', 'main...origin/main') == '0\t0'
    assert git('rev-parse', 'development/retained-dense-step-control-2026-10-06') == BASELINE
    assert git('rev-list', '--left-right', '--count', 'development/retained-dense-step-control-2026-10-06...origin/development/retained-dense-step-control-2026-10-06') == '0\t0'
    assert git('worktree', 'list', '--porcelain').count('worktree ') == 1
    assert read(WORK / 'actual-results-review.json')['accepted']
    assert read(WORK / 'results.json')['accepted']
    ci = read(WORK / 'workflow-final.json')
    assert ci['headSha'] == HEAD and ci['conclusion'] == 'success'
    for case in CASES:
        w = WORK / case
        r = read(w / 'native-control-bridge.json')
        assert r['accepted'] and r['linux_child_birth_absent'] and r['linux_child_group_absent']
        assert birth(r['linux_child_pid']) != r['linux_child_start_ticks']
        assert read(w / 'retirement-run.json')['accepted'] and read(w / 'native-retirement.json')['pass']
        assert read(w / 'observation-analysis.json')['accepted']
    plan = read(WORK / 'execution-plan.json')
    assert plan['source_revision'] == HEAD
    for path, expected in plan['input_seals'].items():
        assert seal(ROOT / path) == expected
    for p in plan['providers']:
        assert seal(Path(p['linux_path'])) == {k: p[k] for k in ('bytes', 'sha256')}
def final_guard():
    assert read(WORK / 'preservation-review.json')['accepted']
    assert read(WORK / 'publication-readback.json')['exact']
    assert read(WORK / 'final-retirement-run.json')['accepted']
    r = read(WORK / 'final-native-retirement.json')
    assert r['pass'] and len(r['targets']) == 6 and not r['matches']
    assert all(t['recorded_birth_absent'] for t in r['targets'])
mode = sys.argv[1]
assert len(sys.argv) == 2 and mode in ('archive', 'close', 'remove')
guard()
if mode == 'archive':
    assert not ARCHIVE.exists()
    ARCHIVE.mkdir(parents=True)
    env = dict(os.environ, GIT_INDEX_FILE=str(WORK / 'preservation.index'))
    useful_names = {'results.json', 'actual-results-review.json', 'observation-analysis.json',
                    'controller-review.json', 'analyzer-source-review.json', 'execution-plan.json',
                    'native-stage-preparation.json', 'source-review-final.json', 'format-token-check.json',
                    'final-retirement-source-review.json'}
    useful = [p for p in files() if p.suffix in ('.py', '.ps1') or p.name in useful_names]
    blobs = {}
    try:
        # Keep the complete reviewed source tree in the tagged evidence snapshot.
        git('read-tree', HEAD, env=env)
        for p in useful:
            name = PREFIX + str(p.relative_to(WORK))
            blob = git('hash-object', '-w', str(p))
            git('update-index', '--add', '--cacheinfo', '100644', blob, name, env=env)
            blobs[str(p.relative_to(WORK))] = {'git_path': name, 'git_blob': blob}
        tree = git('write-tree', env=env)
        message = ('Preserve reviewed native clock stream controls and findings (#44)\n\n'
                   'Diagnostic host/source wiring at 8f1118d preserves original numerical, ordered-stream, refusal and recovery assertions. '
                   'All three original native controls pass; 68 captures independently reproduced. '
                   'Ten Linux authority controls and exact Linux/macOS/Windows non-render CI pass. '
                   'Transport RHS bands dominate these inclusive perturbed raw-clock intervals; no exclusive cost or ordinary benefit is established. '
                   'Scientific rendering, CPU integration and release qualification remain open. '
                   'Full raw/build evidence remains in ignored local attestations, without remote-backup claim.\n')
        commit = git('commit-tree', tree, '-p', HEAD, input=message)
        git('tag', TAG, commit)
    finally:
        for p in WORK.glob('preservation.index*'):
            p.unlink()
    rows = []
    for p in files():
        dest = ARCHIVE / 'files' / p.relative_to(WORK)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, dest)
        assert seal(p) == seal(dest)
        rows.append({'origin': str(p.relative_to(ROOT)), 'recovery': str(dest.relative_to(ROOT)),
                     **seal(p), **blobs.get(str(p.relative_to(WORK)), {})})
    manifest = {'source_revision': HEAD, 'baseline_revision': BASELINE, 'source_tag': TAG,
                'source_commit': commit, 'files': len(rows), 'bytes': sum(r['bytes'] for r in rows),
                'records': rows, 'scope': 'Complete finite WORK recovery; complete diagnostic source plus useful tooling and concise findings in tagged Git. Raw/build evidence remains local. Scoped recorded ownership only; no hidden/global/scheduled absence claim.'}
    (ARCHIVE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({k: manifest[k] for k in ('source_tag', 'source_commit', 'files', 'bytes')}))
elif mode == 'close':
    final_guard()
    manifest = read(ARCHIVE / 'manifest.json')
    assert git('ls-remote', 'origin', 'refs/tags/' + TAG).split()[0] == manifest['source_commit']
    rows = []
    for p in files():
        dest = ARCHIVE / 'files' / p.relative_to(WORK)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            assert seal(dest) == seal(p)
        else:
            shutil.copyfile(p, dest)
        assert seal(dest) == seal(p)
        rows.append({'origin': str(p.relative_to(ROOT)), 'recovery': str(dest.relative_to(ROOT)), **seal(p)})
    for r in manifest['records']:
        assert seal(ROOT / r['recovery']) == {k: r[k] for k in ('bytes', 'sha256')}
    pub = ARCHIVE / 'publication'
    pub.mkdir()
    ledger = {'manifest': seal(ARCHIVE / 'manifest.json'), 'source_tag': TAG,
              'source_commit': manifest['source_commit'], 'records': rows}
    (pub / 'closing-ledger.json').write_text(json.dumps(ledger, indent=2) + '\n')
    print(json.dumps({'files': len(rows), 'bytes': sum(r['bytes'] for r in rows), 'ledger': seal(pub / 'closing-ledger.json')}))
else:
    final_guard()
    ledger = read(ARCHIVE / 'publication/closing-ledger.json')
    assert git('ls-remote', 'origin', 'refs/tags/' + TAG).split()[0] == ledger['source_commit']
    expected = {r['origin']: {k: r[k] for k in ('bytes', 'sha256')} for r in ledger['records']}
    assert {str(p.relative_to(ROOT)): seal(p) for p in files()} == expected
    for r in ledger['records']:
        assert seal(ROOT / r['recovery']) == {k: r[k] for k in ('bytes', 'sha256')}
    shutil.rmtree(WORK)
    state = {'source_revision': git('rev-parse', 'HEAD'), 'clean': not git('status', '--porcelain'),
             'tracking': git('rev-list', '--left-right', '--count', 'HEAD...@{u}'),
             'main_tracking': git('rev-list', '--left-right', '--count', 'main...origin/main'),
             'worktrees': git('worktree', 'list'), 'work_absent': not WORK.exists(),
             'out_entries': sorted(p.name for p in (ROOT / 'out').iterdir()),
             'removed_files': len(expected), 'removed_bytes': sum(v['bytes'] for v in expected.values()),
             'scope': 'Finite diagnostic batch closed; long-term rendering and qualification remediation remains active.'}
    assert state['clean'] and state['tracking'] == state['main_tracking'] == '0\t0'
    (ARCHIVE / 'publication/removal-state.json').write_text(json.dumps(state, indent=2) + '\n')
    print(json.dumps(state))
