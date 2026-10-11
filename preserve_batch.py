"""Preserve useful sources in Git and exact finite WORK artifacts locally."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;A=R/'attestations/native-vulkan/2053ef7/native-shader-clock-admission'
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*args,env=None,input=None):return subprocess.check_output(['git',*args],env=env,input=input,text=True,timeout=30).strip()
assert git('rev-parse','HEAD')=='2053ef79c90752d37b0075b8137cfb9b571b6820'==git('rev-parse','@{u}') and not git('status','--porcelain')
assert not A.exists();A.mkdir(parents=True)
files=sorted(p for p in W.rglob('*') if p.is_file());assert all(not p.is_symlink() for p in files)
# Never copy the temporary Git index or descendants of Windows compiler temp dirs to Git.
useful=[p for p in files if not ({'owner-temp','native-temp','preparation-v0','bridge-pre-review','phase-preparation-v0','failed-phase-compiler-v0'} & set(p.relative_to(W).parts)) and p.suffix in {'.py','.ps1','.cs','.slang'}]
useful += [W/n for n in ('findings.json','source-review.json','source-remediation-review.json','phase-source-review.json','source-validation.json','regression-before-after.json','plan.json','phase-preparation.json')]
assert len(useful)==len(set(useful))
env=dict(os.environ);env['GIT_INDEX_FILE']=str(W/'snapshot.index')
try:
 git('read-tree','--empty',env=env)
 for p in useful:
  blob=git('hash-object','-w',str(p));git('update-index','--add','--cacheinfo','100644',blob,str(p.relative_to(W)),env=env)
 tree=git('write-tree',env=env)
 commit=git('commit-tree',tree,'-p',git('rev-parse','HEAD'),input='Preserve native shader-clock admission sources and finite findings\n\nActual read-only feature support and validated native-default observer modules; no logical-device or GPU clock execution qualification. Production builder correction is2053ef7.\n')
 tag='evidence/issue-44-native-shader-clock-admission-2053ef7'
 git('tag',tag,commit)
finally:
 for p in W.glob('snapshot.index*'):p.unlink()
records=[]
for p in files:
 rel=p.relative_to(W);dest=A/'files'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest);assert seal(p)==seal(dest)
 record={'origin':str(p.relative_to(R)),'recovery':str(dest.relative_to(R)),**seal(p)}
 if p in useful:record['git_recovery']={'commit':commit,'path':str(rel),'blob':git('rev-parse',commit+':'+str(rel))}
 records.append(record)
manifest={'source_revision':git('rev-parse','HEAD'),'native_query_revision':'6be3132b7023f9d6497ab8f67dc779cf8eb4a80c','source_commit':commit,'source_tag':tag,'records':records,'useful_source_files':len(useful),'scope':'Full-byte finiteWORK recovery. Detailed local evidence has no remote-backup claim. No scientific/render/clock-execution/release qualification.'}
(A/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'archive':str(A.relative_to(R)),'manifest':seal(A/'manifest.json'),'files':len(records),'bytes':sum(p['bytes'] for p in records),'source_commit':commit,'source_tag':tag,'useful_source_files':len(useful)}
(W/'preservation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
