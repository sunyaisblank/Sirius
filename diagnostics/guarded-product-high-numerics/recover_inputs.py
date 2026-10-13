"""Recover declared non-Git frozen input versions before retiring this batch."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
revision = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == revision
assert not subprocess.check_output(['git', 'status', '--porcelain'])
tracked = set(subprocess.check_output(['git', 'ls-files', '-z']).decode().split('\0'))


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


records = {}
owners = sorted(work.glob('*/owner.json'))
assert len(owners) == 17
for path in owners:
    owner = json.loads(path.read_text())
    assert owner['status'] == 'terminal' and owner['passed'] and not owner['remaining_owned_processes']
    assert not owner['cleanup_errors'] and all(item['absent'] for item in owner['observed_births'])
    for name, expected in json.loads((path.parent / 'inputs-before.json').read_text()).items():
        if name not in tracked:
            records.setdefault((name, expected['sha256']), expected)
versions = work / 'input-versions'
versions.mkdir(exist_ok=True)
recovery = []
for (name, digest), expected in sorted(records.items()):
    content = {key: expected[key] for key in ['bytes', 'sha256']}
    source = Path(name) if Path(name).is_absolute() else root / name
    target = versions / digest
    if target.exists():
        assert identity(target) == content
        method = 'saved-content-addressed-input'
    else:
        assert identity(source) == content, name
        shutil.copy2(source, target)
        assert identity(source) == content and identity(target) == content
        method = 'unchanged-frozen-input-copy'
    recovery.append({'input': name, **expected, 'saved_path': str(target.relative_to(root)), 'method': method})
report = {'source_revision': revision, 'terminal_owner_records': len(owners),
    'scope': 'Recorded declared per-action frozen input recovery; no claim to freeze every implicit compiler or operating-system dependency.',
    'tracked_source_recovery': 'Original tracked files recoverable from exact unchanged Git source revision; per-action whole declared input hashes retained.',
    'records': recovery}
(work / 'input-recovery.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'non_git_input_bindings_recovered': len(recovery),
    'unique_saved_content_payloads': len(list(versions.iterdir())),
    'unique_saved_bytes': sum(p.stat().st_size for p in versions.iterdir())}))
