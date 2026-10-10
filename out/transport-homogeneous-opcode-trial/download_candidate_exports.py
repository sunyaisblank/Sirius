from pathlib import Path,PurePosixPath
import hashlib,json,stat,subprocess,time,zipfile
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
head=json.loads((WORK/'candidate-source.json').read_text())['source_revision'];run=38048117125
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head and not subprocess.check_output(['git','status','--porcelain'])
meta=json.loads(subprocess.check_output(['gh','api',f'repos/sunyaisblank/Sirius/actions/runs/{run}/artifacts'],text=True,timeout=60))
(WORK/'candidate-native-ci-artifacts.json').write_text(json.dumps(meta,indent=2)+'\n')
for platform in ('windows','macos'):
 matches=[a for a in meta['artifacts'] if a['name']==f'sirius-{platform}-native-attestation'];assert len(matches)==1
 a=matches[0];assert not a['expired'] and a['workflow_run']['id']==run and a['workflow_run']['head_sha']==head
 archive=WORK/f'candidate-{platform}-export.zip';destination=ROOT/'attestations/native-build'/head[:7]/f'{platform}-build'
 assert not archive.exists() and not destination.exists()
 args=['gh','api',f"repos/sunyaisblank/Sirius/actions/artifacts/{a['id']}/zip"];started=time.monotonic()
 with archive.open('wb') as out,(WORK/f'candidate-{platform}-export-download.stderr').open('wb') as err:
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
 (WORK/f'candidate-{platform}-export-download.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ('platform','source_revision','download_bytes','sha256','elapsed_seconds')}),flush=True)
