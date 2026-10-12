"""Authenticate one official diagnostic Windows export; safe fresh extraction."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path,PurePosixPath
import hashlib,json,stat,subprocess,time,zipfile
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
source=json.loads((WORK/'source.json').read_text());head=source['source_revision']
run=json.loads((WORK/'ci-run.json').read_text())['run_id']
ci=json.loads((WORK/'ci.json').read_text());assert ci['accepted'] and ci['source_revision']==head and ci['run_id']==run
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
base=f'repos/sunyaisblank/Sirius/actions/runs/{run}'
meta=json.loads(subprocess.check_output(['gh','api',base+'/artifacts'],timeout=60))
(WORK/'ci-artifacts.json').write_text(json.dumps(meta,indent=2)+'\n')
matches=[a for a in meta['artifacts'] if a['name']=='sirius-windows-native-attestation'];assert len(matches)==1
a=matches[0];assert not a['expired'] and a['workflow_run']['id']==run and a['workflow_run']['head_sha']==head
archive=WORK/'artifact.zip';destination=WORK/'artifact'
assert not archive.exists() and not destination.exists()
argv=['gh','api',f"repos/sunyaisblank/Sirius/actions/artifacts/{a['id']}/zip"]
started=time.monotonic()
with archive.open('xb') as out,(WORK/'artifact-download.stderr').open('xb') as err:
    r=subprocess.run(argv,stdout=out,stderr=err,timeout=180)
assert r.returncode==0
digest=hashlib.sha256(archive.read_bytes()).hexdigest()
assert a['digest']=='sha256:'+digest and archive.stat().st_size==a['size_in_bytes']
with zipfile.ZipFile(archive) as z:
    entries=z.infolist();names=set()
    for e in entries:
        p=PurePosixPath(e.filename)
        assert not p.is_absolute() and '..' not in p.parts and '\\' not in e.filename and ':' not in e.filename
        assert e.filename not in names;names.add(e.filename)
        assert not stat.S_ISLNK(e.external_attr>>16)
    assert len([e for e in entries if not e.is_dir()])==40
    destination.mkdir();z.extractall(destination)
files=[]
for p in sorted(destination.rglob('*')):
    assert not p.is_symlink()
    if p.is_file():files.append({'path':str(p.relative_to(destination)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert len(files)==40
document=json.loads((destination/'windows-build.json').read_text())
assert document['status']=='pass' and document['source_revision']==head and document['domains']==['windows-native-build']
assert document['claims']['runtime_ready'] is False
record={'accepted':True,'source_revision':head,'run_id':run,'artifact':a,'arguments':argv,'returncode':r.returncode,'elapsed_seconds':time.monotonic()-started,'download_bytes':archive.stat().st_size,'sha256':digest,'extracted_files':files,'scope':'Exact official diagnostic build export only. Actual native/scientific execution and full qualification remain pending.'}
(WORK/'artifact-extraction.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('accepted','source_revision','run_id','download_bytes','sha256','elapsed_seconds')}))
