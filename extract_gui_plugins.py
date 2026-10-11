"""Add original packaged Qt plugins before the first profiler help launch."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path, PurePosixPath
import hashlib, json, stat, subprocess, zipfile
root = Path.cwd(); work = Path(__file__).resolve().parent
record = json.loads((work / 'tools-download.json').read_text())
archive = root / record['archive']; target = archive.parent / 'tools'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == record['source_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert archive.stat().st_size == record['bytes'] and hashlib.sha256(archive.read_bytes()).hexdigest() == record['sha256']
assert not (work / 'tools-gui-extraction.json').exists()
files = []
with zipfile.ZipFile(archive) as z:
    for v in z.infolist():
        p = PurePosixPath(v.filename)
        if len(p.parts) < 3 or p.parts[1] != 'plugins' or v.is_dir(): continue
        assert p.parts[0] == 'RadeonDeveloperToolSuite-2026-10-02-1872'
        assert not p.is_absolute() and '..' not in p.parts and '\\' not in v.filename
        assert not stat.S_ISLNK(v.external_attr >> 16)
        q = target.joinpath(*p.parts[1:]); q.parent.mkdir(parents=True, exist_ok=True)
        data = z.read(v); assert len(data) == v.file_size
        with q.open('xb') as f: f.write(data)
        assert q.read_bytes() == data
        files.append({'path': str(q.relative_to(root)), 'zip_member': v.filename, 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(), 'CRC32': v.CRC})
assert any(v['path'].endswith('/plugins/platforms/qwindows.dll') for v in files)
result = {'source_revision': record['source_revision'], 'archive_sha256': record['sha256'], 'files': files, 'file_count': len(files), 'bytes': sum(v['bytes'] for v in files), 'executed_tool_contents': False, 'scope': 'Add exact packaged Qt platform/style/image plugins before the first RGP help launch; no installation, service or executable launch.'}
(work / 'tools-gui-extraction.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ('file_count', 'bytes', 'executed_tool_contents')}))
