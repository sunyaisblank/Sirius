"""Keep publication receipts and readable private patches with the prototype checkpoint."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import stat
import subprocess

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
PARENT = '7e6e97a3e1305c88650cf266f6350a33b534474d'
BRANCH = 'development/software-single-invocation-2026-10-20'
TARGET = ROOT / 'attestations/qualification/630844e7/software-single-invocation-publication'
INDEX = HERE / 'preservation.index'


def git(*args, data=None, env=None):
    return subprocess.check_output(['git', *args], input=data, env=env)


def identity(path):
    info = path.lstat()
    assert stat.S_ISREG(info.st_mode) and not path.is_symlink()
    return {'bytes': info.st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest(), 'mode': oct(stat.S_IMODE(info.st_mode))}


def inventory():
    return {str(path.relative_to(ROOT)): identity(path) for path in sorted(HERE.iterdir()) if path != INDEX}


publication = json.loads((HERE / 'publication.json').read_text())
assert publication['passed'] and publication['prototype_commit'] == PARENT
assert publication['issue_body_and_acceptance_unchanged'] and publication['issue_remains_open']
assert git('rev-parse', 'HEAD').decode().strip() == publication['source_revision'] and git('status', '--porcelain') == b''
assert git('rev-parse', 'refs/heads/' + BRANCH).decode().strip() == PARENT
assert not TARGET.exists() and not INDEX.exists()
before = inventory()
TARGET.mkdir(parents=True)
copies = {}
for name, row in before.items():
    destination = TARGET / 'files' / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / name, destination)
    assert identity(destination) == row
    copies['files/' + name] = {'source': name, **row}
for name, source in {
    'prototype.patch': 'out/software-scalar-nojoin/candidate.patch',
    'producer-include.patch': 'out/software-scalar-numerics-build/source.patch',
    'prototype-source-contracts.json': 'out/software-scalar-nojoin/source-contracts.json',
}.items():
    path = ROOT / source
    shutil.copy2(path, TARGET / name)
    assert identity(TARGET / name) == identity(path)
    copies[name] = {'source': source, **identity(path)}
manifest = {'schema': 'sirius-software-scalar-prototype-publication-v1', 'operator_source_revision': publication['source_revision'],
            'prototype_first_checkpoint': PARENT, 'prototype_branch': BRANCH, 'payloads': copies,
            'issue_comment': publication['comment_url'], 'issue44_acceptance_unchanged': True,
            'source_working_tree_clean_one_checkout_no_live_owned_processes': True,
            'unfinished_private_prototype': True, 'prototype_adopted': False,
            'source_inventory_before_after_equal': True,
            'adoption_boundary': 'The isolated producer was forced to Lavapipe. Its replacement of one normal-array is not a production routing change. Any future production integration must preserve all existing nonsoftware modules/routes, exact precision admission and original controls/benefit gates; no result here admits MoltenVK/Dozen/other physical drivers.',
            'workspace': 'The original283-file prototype snapshot and two private bin stages remain immutable inputs for continuing work. Only this finished publication scratch is disposable after exact Git recovery; preserve canonical bins and unrelated files.'}
(TARGET / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
assert before == inventory()
environment = dict(os.environ, GIT_INDEX_FILE=str(INDEX))
try:
    git('read-tree', PARENT, env=environment)
    rows = []
    for path in sorted(TARGET.rglob('*')):
        if path.is_file():
            blob = git('hash-object', '-w', str(path)).decode().strip()
            rows.append(('100644 ' + blob + '\t' + str(path.relative_to(ROOT)) + '\n').encode())
    git('update-index', '--index-info', data=b''.join(rows), env=environment)
    tree = git('write-tree', env=environment).decode().strip()
    commit = git('commit-tree', tree, '-p', PARENT, data=b'Publish scalar Endpoint prerequisite evidence and preserve readable prototype patches\n').decode().strip()
    git('update-ref', 'refs/heads/' + BRANCH, commit, PARENT)
    for path in TARGET.rglob('*'):
        if path.is_file():
            data = git('show', commit + ':' + str(path.relative_to(ROOT)))
            row = identity(path)
            assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
finally:
    INDEX.unlink(missing_ok=True)
assert before == inventory() and git('status', '--porcelain') == b''
print(json.dumps({'commit': commit, 'branch': BRANCH, 'bundle': str(TARGET.relative_to(ROOT)), 'payloads': len(copies)}))
