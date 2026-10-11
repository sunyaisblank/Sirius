"""Acquire the pinned official distribution without executing its contents."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, subprocess, time, urllib.request, zipfile
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
TARGET = ROOT / 'bin/windows-msvc/qualification-tools/rdts-2026-10-02-1872'
URL = 'https://gpuopen.com/download/RadeonDeveloperToolSuite-2026-10-02-1872.zip'
EXPECTED_BYTES = 424481206
def git(*a): return subprocess.check_output(['git', *a], text=True).strip()
head = json.loads((WORK / 'source.json').read_text())['source_revision']
assert git('rev-parse','HEAD') == head == git('rev-parse','@{u}') and not git('status','--porcelain')
assert not TARGET.exists()
partial = WORK / 'rdts-download.partial'; assert not partial.exists()
start = time.monotonic(); digest = hashlib.sha256(); count = 0
with urllib.request.urlopen(URL, timeout=15) as r, partial.open('xb') as f:
    assert r.status == 200 and r.url == URL and int(r.headers['Content-Length']) == EXPECTED_BYTES
    headers = dict(r.headers)
    while True:
        assert time.monotonic() - start < 180, 'Official tools download deadline'
        b = r.read(1024 * 1024)
        if not b: break
        count += len(b); assert count <= EXPECTED_BYTES
        digest.update(b); f.write(b)
assert count == EXPECTED_BYTES
with zipfile.ZipFile(partial) as z:
    inventory = [{'name':v.filename,'bytes':v.file_size,'compressed_bytes':v.compress_size,'CRC32':v.CRC} for v in z.infolist()]
assert git('rev-parse','HEAD') == head == git('rev-parse','@{u}') and not git('status','--porcelain')
TARGET.mkdir(parents=True)
archive = TARGET / 'RadeonDeveloperToolSuite-2026-10-02-1872.zip'
os.rename(partial, archive)
record = {'source_revision':head,'url':URL,'resolved_url':URL,'headers':headers,'archive':str(archive.relative_to(ROOT)),'bytes':count,'sha256':digest.hexdigest(),'wall_seconds':time.monotonic()-start,'inventory':inventory,'executed_tool_contents':False,'meaning':'Official HTTPS archive identity and inventory only; not a native profiling or qualification result.'}
(WORK/'tools-download.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('archive','bytes','sha256','wall_seconds','executed_tool_contents')}))
print(json.dumps(sorted(set(v['name'].split('/')[0] for v in inventory))))
