import hashlib
import json
import os
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PARENT = '506fc0622a03b57c7c5965421e2bb48120cbbc79'
BRANCH = 'refs/heads/development/native-fp32-packed-scratch-runtime-2026-10-18'
TAG = 'refs/tags/evidence/native-coefficient-sharing-compiler-driver-2c8ef086'
BUNDLE = ROOT / 'attestations/qualification/2c8ef086/native-coefficient-sharing-source-compiler-driver'


def run(*args, data=None, env=None):
    return subprocess.run(args, cwd=ROOT, input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True, env=env).stdout


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'mode': oct(stat.S_IMODE(path.stat().st_mode))}


assert not run('git', 'status', '--porcelain')
assert run('git', 'rev-parse', BRANCH).decode().strip() == PARENT
assert not BUNDLE.exists()
for name in ['conservation.json', 'retirement.json', 'published-input-retirement.json']:
    assert json.loads((ROOT / 'out/native-rounded-error-closeout' / name).read_text())['passed']
assert json.loads((ROOT / 'out/native-coefficient-cost/driver-cost.json').read_text())['passed']
assert json.loads((ROOT / 'out/native-coefficient-sharing-emitted/emitted-proof.json').read_text())['passed']
assert json.loads((ROOT / 'out/native-coefficient-sharing-compiler/candidate/owner.json').read_text())['passed']
assert json.loads((WORK / 'source-Git.json').read_text())['fresh_remote_verified']
roots = ['out/native-coefficient-sharing-' + name for name in
         ['source', 'compiler', 'emitted', 'driver', 'producer-source', 'runtime-source', 'ci-publication']]
roots += ['out/native-coefficient-cost', 'out/native-rounded-error-closeout', 'out/native-rounded-error-preservation']
payloads = {}
for name in roots:
    for source in sorted((ROOT / name).rglob('*')):
        assert not source.is_symlink()
        if not source.is_file():
            continue
        assert source.suffix != '.index'
        rel = str(source.relative_to(ROOT))
        target = BUNDLE / 'files' / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        before = pin(source)
        shutil.copyfile(source, target)
        target.chmod(stat.S_IMODE(source.stat().st_mode))
        assert before == pin(source) == pin(target)
        payloads['files/' + rel] = {'source': rel, **before}
manifest = {'schema': 'sirius-native-coefficient-sharing-compiler-driver-preservation-v1',
            'canonical_revision': '2c8ef08693d46f57d5dea14efdbddd4874dac684',
            'private_producer_revision': 'e1c618f83861c478fe18c5a22f5d9aa27f6280fd',
            'payloads': payloads,
            'milestones': ['Complete rejected rounded-error conservation, unchanged issue44 acceptance and exact owned scratch retirement',
                           'Distinct coefficient-call selection source proof, one original Transport compiler screen and exact emitted clone/ABI proof',
                           'One zero-dispatch actual AMD driver query with complete ISA/text joins: text22.0701% smaller, VGPR55to63 countercost',
                           'Exact private source published and sole compile-only Windows CI in progress; original runtime source held false/artifact pending'],
            'runtime_admitted': False, 'source_adopted': False, 'full_qualification_claimed': False,
            'external_sources': [{'commit': PARENT, 'scope': 'Complete rounded-error rejected raw runtime/ZIP and prior source/compiler ancestry'},
                                 {'commit': 'e1c618f83861c478fe18c5a22f5d9aa27f6280fd', 'scope': 'Complete actual77 producer payloads/composition and Windows apparatus/source; not copied again'},
                                 {'commit': '1cb6e6aaff04ddd54ff7ee133b773257338f2c52', 'scope': 'Saved original Transport driver object and exact decoder; baseline not re-executed'}],
            'limits': 'No new provider arithmetic, eight-control, benefit, frame, installed-release or science verdict. Diagnostic Vulkan1.2 pipeline differs from application Vulkan1.3 admission. Host input seals are bounded hash records, not backups.'}
(BUNDLE / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
index = WORK / 'milestone.index'
assert not index.exists()
env = dict(os.environ, GIT_INDEX_FILE=str(index))
records = {}
try:
    run('git', 'read-tree', PARENT, env=env)
    rows = []
    for path in sorted(BUNDLE.rglob('*')):
        if not path.is_file():
            continue
        data = path.read_bytes()
        blob = run('git', 'hash-object', '-w', '--stdin', data=data).decode().strip()
        name = str(path.relative_to(ROOT))
        mode = '100755' if path.stat().st_mode & 0o111 else '100644'
        records[name] = {**pin(path), 'Git_mode': mode}
        rows.append(f'{mode} {blob}\t{name}\0'.encode())
    run('git', 'update-index', '-z', '--index-info', data=b''.join(rows), env=env)
    tree = run('git', 'write-tree', env=env).decode().strip()
    commit = run('git', 'commit-tree', tree, '-p', PARENT,
                 data=b'Conserve rounded-error rejection and coefficient-sharing compiler/driver reduction\n').decode().strip()
    assert set(run('git', 'diff-tree', '--no-commit-id', '--name-only', '-r', commit).decode().splitlines()) == set(records)
    for name, row in records.items():
        data = run('git', 'show', commit + ':' + name)
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
        assert run('git', 'ls-tree', commit, '--', name).decode().split()[0] == row['Git_mode']
    run('git', 'update-ref', BRANCH, commit, PARENT)
    run('git', 'update-ref', TAG, commit, '0' * 40)
    receipt = {'schema': 'sirius-native-coefficient-sharing-milestone-Git-v1',
               'commit': commit, 'parent': PARENT, 'branch': BRANCH, 'tag': TAG,
               'bundle': str(BUNDLE.relative_to(ROOT)), 'files': records,
               'payload_files': len(payloads), 'payload_bytes': sum(row['bytes'] for row in payloads.values()),
               'stored_full_bytes_and_Git_modes_verified': True, 'pushed': False,
               'runtime_admitted': False, 'source_adopted': False}
    (WORK / 'milestone-Git.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'files'}))
finally:
    index.unlink(missing_ok=True)
