import hashlib
import json
import os
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PARENT = 'a2c34d8db649533ce8bf853ebeeaaada35c81a07'
BRANCH = 'refs/heads/development/native-fp32-packed-scratch-runtime-2026-10-18'
TAG = 'refs/tags/evidence/native-rounded-error-source-compiler-2c8ef086'
BUNDLE = ROOT / 'attestations/qualification/2c8ef086/native-rounded-error-source-compiler'


def run(*args, data=None, env=None):
    return subprocess.run(args, cwd=ROOT, input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True, env=env).stdout


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'mode': oct(stat.S_IMODE(path.stat().st_mode))}


assert run('git', 'rev-parse', 'HEAD').decode().strip() == '2c8ef08693d46f57d5dea14efdbddd4874dac684'
assert not run('git', 'status', '--porcelain')
assert run('git', 'rev-parse', BRANCH).decode().strip() == PARENT
assert not BUNDLE.exists()
for name in ['root-source-admission.json', 'root-compiler-terminal.json',
             'root-emitted-admission.json', 'root-producer-source-admission.json']:
    record = json.loads((ROOT / 'out/native-rounded-error-publication' / name).read_text())
    assert record.get('passed', record.get('admitted_for_narrow_compiler_screen')) is True
payloads = {}
for prefix in ['out/native-rounded-error-source', 'out/native-rounded-error-compiler-source',
               'out/native-rounded-error-emitted-review', 'out/native-rounded-error-producer-source',
               'out/native-rounded-error-pe-audit', 'out/native-rounded-error-runtime-source',
               'out/native-rounded-error-publication']:
    for source in sorted((ROOT / prefix).rglob('*')):
        assert not source.is_symlink()
        if not source.is_file():
            continue
        name = str(source.relative_to(ROOT))
        target = BUNDLE / 'files' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        target.chmod(stat.S_IMODE(source.stat().st_mode))
        assert pin(source) == pin(target)
        payloads['files/' + name] = {'source': name, **pin(source)}
source = Path(__file__).resolve()
target = BUNDLE / 'files' / str(source.relative_to(ROOT))
target.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(source, target)
target.chmod(stat.S_IMODE(source.stat().st_mode))
payloads[str(target.relative_to(BUNDLE))] = {'source': str(source.relative_to(ROOT)), **pin(source)}
manifest = {'schema': 'sirius-native-rounded-error-source-compiler-preservation-v1',
            'canonical_source_revision': '2c8ef08693d46f57d5dea14efdbddd4874dac684',
            'private_producer_source_revision': '625abe4569c4b377fad6d8586d6f0487f01d9903',
            'payloads': payloads,
            'scope': 'Exact reviewed bounded native error shortcut source, one successful five-module compiler '
                     'screen, full emitted comparison and held producer/runtime source. No actual Windows '
                     'producer, provider correctness, benefit, adoption or full qualification verdict.',
            'external_preserved_references': [
                {'commit': '625abe4569c4b377fad6d8586d6f0487f01d9903',
                 'scope': 'Published complete private95 additions including five candidate source files, '
                          '77 payloads/composition, nine Windows source files, root emitted admission and '
                          'source publisher/workflow; all538 canonical files remain unchanged.'},
                {'commit': 'fd530d957d2e0737f5453d88d78d979f8ed17607',
                 'scope': 'Original source/compiler/runtime preimages and adopted constructor baseline.'},
                {'commit': '1cb6e6aaff04ddd54ff7ee133b773257338f2c52',
                 'scope': 'Saved default native driver cost premise only; wave32 variant remains rejected.'}],
            'limits': 'Compiler/tool/host dependency seals are bounded hash records, not backups of '
                      'external tools/providers. Runtime protocol IDs are deliberately pending in this '
                      'source snapshot; later data finalization requires its own record.'}
(BUNDLE / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
index = WORK / 'source-compiler.index'
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
        mode = '100755' if path.stat().st_mode & 0o111 else '100644'
        name = str(path.relative_to(ROOT))
        records[name] = {**pin(path), 'Git_mode': mode}
        rows.append(f'{mode} {blob}\t{name}\0'.encode())
    run('git', 'update-index', '-z', '--index-info', data=b''.join(rows), env=env)
    tree = run('git', 'write-tree', env=env).decode().strip()
    commit = run('git', 'commit-tree', tree, '-p', PARENT,
                 data=b'Preserve rounded-product error exact source and emitted compiler proof\n').decode().strip()
    for name, row in records.items():
        data = run('git', 'show', commit + ':' + name)
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
    run('git', 'update-ref', BRANCH, commit, PARENT)
    run('git', 'update-ref', TAG, commit, '0' * 40)
    receipt = {'schema': 'sirius-native-rounded-error-source-compiler-Git-v1', 'commit': commit,
               'parent': PARENT, 'branch': BRANCH, 'tag': TAG,
               'bundle': str(BUNDLE.relative_to(ROOT)), 'files': records,
               'full_Git_bytes_verified': True, 'pushed': False,
               'runtime_admitted': False, 'source_adopted': False}
    (WORK / 'source-compiler-Git.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'files'}))
    print(json.dumps({'payload_files': len(payloads), 'Git_files': len(records),
                      'payload_bytes': sum(row['bytes'] for row in payloads.values())}))
finally:
    index.unlink(missing_ok=True)
