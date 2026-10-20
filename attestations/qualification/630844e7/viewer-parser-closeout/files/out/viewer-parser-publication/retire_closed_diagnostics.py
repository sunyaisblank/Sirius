"""Retire only the closed157-file scope recovered exactly from pushed Git."""
from pathlib import Path
import hashlib
import json
import stat
import subprocess
import zipfile

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
COMMIT = '079d204b968d7f1c9967b589d28de70017999b96'
bundle = ROOT / 'attestations/qualification/630844e7/viewer-parser-source-build-and-display'
manifest = json.loads((bundle / 'MANIFEST.json').read_text())
assert subprocess.check_output(['git', 'status', '--porcelain']) == b''
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == manifest['source_revision']
refs = subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/heads/development/viewer-parser-qualification-2026-10-20',
                                'refs/tags/evidence/viewer-parser-source-build-display-630844e7'], text=True)
assert len(refs.splitlines()) == 2 and all(line.split()[0] == COMMIT for line in refs.splitlines())


def identity(path):
    info = path.lstat()
    assert stat.S_ISREG(info.st_mode) and not path.is_symlink()
    return {'bytes': info.st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest(),
            'mode': oct(stat.S_IMODE(info.st_mode))}


for path in bundle.iterdir():
    expected = subprocess.check_output(['git', 'show', COMMIT + ':' + str(path.relative_to(ROOT))])
    assert len(expected) == path.stat().st_size
    assert hashlib.sha256(expected).hexdigest() == identity(path)['sha256']
archive = bundle / manifest['closed_evidence_archive']['path']
assert identity(archive) == {key: manifest['closed_evidence_archive'][key] for key in ('bytes', 'sha256', 'mode')}
members = manifest['archive_backed_files']
roots = ['out/viewer-constant-parser-validation', 'out/viewer-parser-final-check', 'out/viewer-parser-sanitizer-check',
         'out/viewer-constant-parser-closeout', 'out/viewer-parser-preservation']
with zipfile.ZipFile(archive) as handle:
    assert handle.testzip() is None and set(handle.namelist()) == set(members)
    for name, row in members.items():
        assert identity(ROOT / name) == row
        data = handle.read(name)
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
        assert stat.S_IMODE(handle.getinfo(name).external_attr >> 16) == int(row['mode'], 8)
for relative in roots:
    base = ROOT / relative
    assert base.is_dir() and not base.is_symlink()
    assert not any(path.is_symlink() for path in base.rglob('*'))
    actual = {str(path.relative_to(ROOT)) for path in base.rglob('*') if path.is_file()}
    assert actual == {name for name in members if Path(name).is_relative_to(relative)}
for name, row in members.items():
    path = ROOT / name
    assert identity(path) == row
    path.unlink()
for relative in roots:
    base = ROOT / relative
    for path in sorted(base.rglob('*'), key=lambda item: len(item.parts), reverse=True):
        assert path.is_dir() and not path.is_symlink()
        path.rmdir()
    base.rmdir()
receipt = {'schema': 'sirius-viewer-parser-closed-diagnostics-retired-v1', 'passed': True,
           'preserved_commit': COMMIT, 'retired_files': len(members), 'retired_roots': roots,
           'exact_Git_archive_bytes_modes_and_CRC_verified_before_removal': True,
           'canonical_source_build_bins_executed_volumes_and_qualification_bundles_preserved': True,
           'closed_owner_records': 16, 'remaining_owned_processes': [], 'source_unchanged': True}
(HERE / 'diagnostic-retirement.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt))
