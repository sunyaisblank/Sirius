"""Publish reviewed bodies only if the remote originals still match."""
from pathlib import Path
import json
import subprocess

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent


def query(kind, number):
    fields = 'body,state,updatedAt,url' + (',headRefOid,headRefName,baseRefName,isDraft' if kind == 'pr' else '')
    return json.loads(subprocess.check_output(['gh', kind, 'view', str(number), '--json', fields]))


prepared = json.loads((HERE / 'prepared.json').read_text())
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == prepared['source_revision']
assert subprocess.check_output(['git', 'status', '--porcelain']) == b''
remote = subprocess.check_output(['git', 'ls-remote', 'origin',
                                 'refs/heads/development/sanitizer-shader-constant-parser-2026-10-20',
                                 'refs/heads/development/viewer-parser-qualification-2026-10-20',
                                 'refs/tags/evidence/viewer-parser-source-build-display-630844e7',
                                 'refs/heads/main'], text=True)
refs = {line.split()[1]: line.split()[0] for line in remote.splitlines()}
assert refs['refs/heads/main'] == '5e706073a3fd2f9afcf8efbcde2f094122050f6a'
assert refs['refs/heads/development/sanitizer-shader-constant-parser-2026-10-20'] == prepared['source_revision']
assert refs['refs/heads/development/viewer-parser-qualification-2026-10-20'] == prepared['evidence_commit']
assert refs['refs/tags/evidence/viewer-parser-source-build-display-630844e7'] == prepared['evidence_commit']
expected = []
for kind, number in [('issue', 44), ('issue', 48), ('issue', 53), ('pr', 62)]:
    original = json.loads((HERE / f'{kind}-{number}-before.json').read_text())
    current = query(kind, number)
    assert current['body'] == original['body'] and current['updatedAt'] == original['updatedAt']
    assert current['state'] == 'OPEN'
    if kind == 'pr':
        assert current['isDraft'] and current['headRefOid'] == '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
        assert current['headRefName'] == 'development/retained-dense-step-control-2026-10-06' and current['baseRefName'] == 'main'
    expected.append((kind, number, original))
    (HERE / f'{kind}-{number}-verified-before.json').write_text(json.dumps(current, indent=2) + '\n')

records = []
for kind, number, original in expected:
    path = HERE / f'{kind}-{number}-body.txt'
    result = subprocess.check_output(['gh', kind, 'edit', str(number), '--body-file', str(path)], text=True).strip()
    actual = query(kind, number)
    assert actual['body'] == path.read_text() and actual['state'] == 'OPEN'
    if number == 44:
        assert actual['body'][actual['body'].index('Acceptance remains:'):] == original['body'][original['body'].index('Acceptance remains:'):]
    if number == 48:
        assert actual['body'][actual['body'].index('Required acceptance:'):] == original['body'][original['body'].index('Required acceptance:'):]
    if kind == 'pr':
        assert all(actual[key] == original[key] for key in ('headRefOid', 'headRefName', 'baseRefName', 'isDraft'))
    (HERE / f'{kind}-{number}-after.json').write_text(json.dumps(actual, indent=2) + '\n')
    records.append({'kind': kind, 'number': number, 'url': actual['url'], 'updated_at': actual['updatedAt'],
                    'expected_body_readback_exact': True, 'state': actual['state']})
    print(kind, number, result)
final = {'schema': 'sirius-viewer-parser-queue-publication-v1', 'passed': True,
         'source_revision': prepared['source_revision'], 'evidence_commit': prepared['evidence_commit'],
         'remote_refs': refs, 'updates': records, 'acceptance_preserved': True,
         'PR_head_base_draft_unchanged': True, 'source_or_main_not_modified': True}
(HERE / 'publication.json').write_text(json.dumps(final, indent=2) + '\n')
