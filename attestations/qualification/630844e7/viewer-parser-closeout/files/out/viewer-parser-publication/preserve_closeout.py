"""Preserve exact queue readbacks and scoped retirement without moving source HEAD."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import stat
import subprocess

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
PARENT = '079d204b968d7f1c9967b589d28de70017999b96'
BRANCH = 'development/viewer-parser-qualification-2026-10-20'
TARGET = ROOT / 'attestations/qualification/630844e7/viewer-parser-closeout'
index = HERE / 'closeout.index'


def git(*args, data=None, env=None):
    return subprocess.check_output(['git', *args], input=data, env=env)


def identity(path):
    info = path.lstat()
    assert stat.S_ISREG(info.st_mode) and not path.is_symlink()
    return {'bytes': info.st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest(),
            'mode': oct(stat.S_IMODE(info.st_mode))}


def inventory():
    return {str(path.relative_to(ROOT)): identity(path) for path in sorted(HERE.iterdir()) if path != index}


publication = json.loads((HERE / 'publication.json').read_text())
retirement = json.loads((HERE / 'diagnostic-retirement.json').read_text())
assert publication['passed'] and retirement['passed'] and retirement['retired_files'] == 157
assert publication['acceptance_preserved'] and publication['PR_head_base_draft_unchanged']
assert git('rev-parse', 'HEAD').decode().strip() == publication['source_revision'] and git('status', '--porcelain') == b''
assert git('rev-parse', 'refs/heads/' + BRANCH).decode().strip() == PARENT
assert not TARGET.exists() and not index.exists()
before = inventory()
TARGET.mkdir(parents=True)
copies = {}
for name, row in before.items():
    destination = TARGET / 'files' / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / name, destination)
    assert identity(destination) == row
    copies['files/' + name] = {'source': name, **row}
manifest = {'schema': 'sirius-viewer-parser-publication-closeout-v1', 'source_revision': publication['source_revision'],
            'source_branch': 'development/sanitizer-shader-constant-parser-2026-10-20',
            'payloads': copies, 'source_inventory_before_after_equal': True,
            'completed_build_display_evidence_commit': PARENT, 'published_issues': [44, 48, 53],
            'all_issues_open_PR62_open_draft_unmerged': True, 'acceptance_unchanged': True,
            'retired_closed_diagnostic_files': 157, 'source_build_checks_and_six_original_display_cases_passed': True,
            'current_exact_source_external_domains_admitted': 0, 'required_domains': 8,
            'goal_active_full_scientific_Mandatory_runtime_and_release_qualification_incomplete': True}
(TARGET / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
assert before == inventory()
environment = dict(os.environ, GIT_INDEX_FILE=str(index))
try:
    git('read-tree', PARENT, env=environment)
    rows = []
    for path in sorted(TARGET.rglob('*')):
        if path.is_file():
            blob = git('hash-object', '-w', str(path)).decode().strip()
            rows.append(('100644 ' + blob + '\t' + str(path.relative_to(ROOT)) + '\n').encode())
    git('update-index', '--index-info', data=b''.join(rows), env=environment)
    tree = git('write-tree', env=environment).decode().strip()
    commit = git('commit-tree', tree, '-p', PARENT, data=b'Record parser qualification issue updates and exact scratch retirement\n').decode().strip()
    git('update-ref', 'refs/heads/' + BRANCH, commit, PARENT)
    for path in TARGET.rglob('*'):
        if path.is_file():
            data = git('show', commit + ':' + str(path.relative_to(ROOT)))
            row = identity(path)
            assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
finally:
    index.unlink(missing_ok=True)
assert before == inventory() and git('status', '--porcelain') == b''
print(json.dumps({'commit': commit, 'branch': BRANCH, 'bundle': str(TARGET.relative_to(ROOT)),
                  'payload_count': len(copies), 'source_checkout_unchanged': True}))
