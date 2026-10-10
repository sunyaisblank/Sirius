from pathlib import Path,PurePosixPath
import hashlib,json,stat,subprocess,time,zipfile
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
label, = __import__('sys').argv[1:]
assert label in ('baseline','candidate')
source=json.loads((WORK/'source.json').read_text())
head=source['baseline_revision'] if label=='baseline' else source['source_revision']
dispatch=json.loads((WORK/(label+'-native-export-dispatch.json')).read_text())
assert dispatch['source_revision']==head and dispatch['returncode']==0
run=int(dispatch['stdout'].strip().rsplit('/',1)[1])
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==source['source_revision'] and not subprocess.check_output(['git','status','--porcelain'])
run_meta=json.loads(subprocess.check_output(['gh','api',f'repos/sunyaisblank/Sirius/actions/runs/{run}'],text=True,timeout=60))
assert run_meta['head_sha']==head and run_meta['status']=='completed' and run_meta['conclusion']=='success'
(WORK/(label+'-native-export-ci-readback.json')).write_text(json.dumps(run_meta,indent=2)+'\n')
meta=json.loads(subprocess.check_output(['gh','api',f'repos/sunyaisblank/Sirius/actions/runs/{run}/artifacts'],text=True,timeout=60))
(WORK/(label+'-native-ci-artifacts.json')).write_text(json.dumps(meta,indent=2)+'\n')
for platform in ('windows','macos'):
 matches=[a for a in meta['artifacts'] if a['name']==f'sirius-{platform}-native-attestation'];assert len(matches)==1
 a=matches[0];assert not a['expired'] and a['workflow_run']['id']==run and a['workflow_run']['head_sha']==head
 archive=WORK/f'{label}-{platform}-export.zip';destination=ROOT/'attestations/native-build'/head[:7]/f'{platform}-build'
 assert not archive.exists() and not destination.exists()
 args=['gh','api',f"repos/sunyaisblank/Sirius/actions/artifacts/{a['id']}/zip"];started=time.monotonic()
 with archive.open('wb') as out,(WORK/f'{label}-{platform}-export-download.stderr').open('wb') as err:
  c=subprocess.run(args,stdout=out,stderr=err,timeout=180)
 assert c.returncode==0
 digest=hashlib.sha256(archive.read_bytes()).hexdigest();assert a['digest']=='sha256:'+digest and archive.stat().st_size==a['size_in_bytes']
 with zipfile.ZipFile(archive) as z:
  entries=z.infolist();names=set()
  for e in entries:
   p=PurePosixPath(e.filename);assert not p.is_absolute() and '..' not in p.parts and '\\' not in e.filename and ':' not in e.filename
   assert e.filename not in names;names.add(e.filename)
   assert not stat.S_ISLNK(e.external_attr>>16)
  assert len([e for e in entries if not e.is_dir()])==40
  destination.mkdir(parents=True);z.extractall(destination)
 files=[]
 for p in sorted(destination.rglob('*')):
  assert not p.is_symlink()
  if p.is_file():files.append({'path':str(p.relative_to(destination)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 assert len(files)==40
 result={'source_revision':head,'platform':platform,'artifact':a,'arguments':args,'returncode':c.returncode,'elapsed_seconds':time.monotonic()-started,'download_bytes':archive.stat().st_size,'sha256':digest,'extracted_files':files,'scope':'Exact official build artifact preservation; no native numerical execution/fullqualification claim.'}
 (WORK/f'{label}-{platform}-export-download.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ('platform','source_revision','download_bytes','sha256','elapsed_seconds')}),flush=True)
