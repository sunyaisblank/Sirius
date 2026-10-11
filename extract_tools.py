"""Extract only profiler/controller dependencies and their exact packaged manuals."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path, PurePosixPath
import hashlib, json, shutil, stat, subprocess, zipfile
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
record = json.loads((WORK / 'tools-download.json').read_text())
archive = ROOT / record['archive']; target = archive.parent / 'tools'
head = json.loads((WORK / 'source.json').read_text())['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip() == head
assert not subprocess.check_output(['git','status','--porcelain'])
assert archive.stat().st_size == record['bytes'] and hashlib.sha256(archive.read_bytes()).hexdigest() == record['sha256']
assert not target.exists()
executables = {'RadeonDeveloperPanelCLI.exe','RadeonDeveloperServiceCLI.exe','RadeonGPUProfiler.exe','rtda.exe'}
files = []
with zipfile.ZipFile(archive) as z:
    assert len(z.namelist()) == len(set(z.namelist()))
    for v in z.infolist():
        p = PurePosixPath(v.filename)
        assert not p.is_absolute() and '..' not in p.parts and '\\' not in v.filename
        assert p.parts[0] == 'RadeonDeveloperToolSuite-2026-10-02-1872'
        assert not stat.S_ISLNK(v.external_attr >> 16)
        relative = PurePosixPath(*p.parts[1:])
        selected = (len(relative.parts) == 1 and (relative.suffix == '.dll' or relative.name in executables or relative.suffix == '.txt')) or (str(relative).startswith(('help/rdp/_sources/','help/rgp/_sources/')) and relative.name.endswith('.rst.txt'))
        if not selected or v.is_dir(): continue
        q = target.joinpath(*relative.parts); q.parent.mkdir(parents=True,exist_ok=True)
        with z.open(v) as inp, q.open('xb') as out: shutil.copyfileobj(inp,out,1024*1024)
        assert q.stat().st_size == v.file_size
        files.append({'path':str(q.relative_to(ROOT)),'zip_member':v.filename,'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'CRC32':v.CRC})
assert executables <= {Path(v['path']).name for v in files}
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip() == head
assert not subprocess.check_output(['git','status','--porcelain'])
output = {'source_revision':head,'official_archive':{k:record[k] for k in ('archive','url','bytes','sha256')},'files':files,'file_count':len(files),'bytes':sum(v['bytes'] for v in files),'executed_tool_contents':False,'scope':'Exact profiler/controller dependencies and original vendor documentation extraction only; no installation/service/driver setting or executable launch.'}
(WORK / 'tools-extraction.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'file_count':len(files),'bytes':output['bytes'],'executed_tool_contents':False}))
