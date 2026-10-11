"""Preserve late closure sources and remove only fully recovered, retired WORK."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;A=R/'attestations/native-vulkan/2053ef7/native-shader-clock-admission';P=A/'publication';P.mkdir(exist_ok=True)
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*args,env=None,input=None):return subprocess.check_output(['git',*args],env=env,input=input,text=True,timeout=30).strip()
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
def guards():
 assert W==R/'out/native-shader-clock-admission' and not W.is_symlink()
 assert git('rev-parse','HEAD')=='2053ef79c90752d37b0075b8137cfb9b571b6820'==git('rev-parse','@{u}') and not git('status','--porcelain')
 assert git('rev-list','--left-right','--count','main...origin/main')=='0\t0'
 assert len(git('worktree','list','--porcelain').split('worktree '))-1==1
 assert read(W/'preservation-review.json')['accepted'] and read(W/'publication-readback.json')['comment']['exact_body_readback']
 assert read(W/'final-retirement-run.json')['accepted'] and read(W/'final-native-retirement.json')['pass']
 for receipt in [W/'compiler-admission.json',W/'phase-compiler.json',W/'failed-phase-compiler-v0/phase-compiler.json']:
  for run in read(receipt)['runs']:
   assert birth(run['pid'])!=run['start_ticks']
 assert birth(read(W/'clock-query-bridge.json')['linux_child_pid'])!=read(W/'clock-query-bridge.json')['linux_child_start_ticks']
 for name,v in read(W/'original-products-before-build.json').items():assert seal(R/name)==v
 for provider in read(W/'plan.json')['providers']:assert seal(Path(provider['linux_path']))=={k:provider[k] for k in ('bytes','sha256')}
 for record in read(A/'manifest.json')['records']:assert seal(R/record['recovery'])=={k:record[k] for k in ('bytes','sha256')}
mode=sys.argv[1];assert mode in ('prepare','remove');guards()
if mode=='prepare':
 env=dict(os.environ);env['GIT_INDEX_FILE']=str(W/'closing.index')
 prior=read(A/'manifest.json')['source_commit']
 try:
  git('read-tree',prior,env=env)
  for n in ('finish_batch.py','run_final_native_retirement.py','preservation-review.json'):
   blob=git('hash-object','-w',str(W/n));git('update-index','--add','--cacheinfo','100644',blob,n,env=env)
  tree=git('write-tree',env=env);commit=git('commit-tree',tree,'-p',prior,input='Preserve final native clock admission recovery and scoped retirement helpers\n')
  tag='evidence/issue-44-native-shader-clock-admission-2053ef7-closed';git('tag',tag,commit)
 finally:
  for p in W.glob('closing.index*'):p.unlink()
 rows=[]
 for p in sorted(p for p in W.rglob('*') if p.is_file()):
  assert not p.is_symlink();dest=A/'files'/p.relative_to(W);dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():assert seal(dest)==seal(p)
  else:shutil.copyfile(p,dest)
  rows.append({'origin':str(p.relative_to(R)),'recovery':str(dest.relative_to(R)),**seal(p)})
 (P/'closing-ledger.json').write_text(json.dumps({'source_revision':git('rev-parse','HEAD'),'initial_manifest':seal(A/'manifest.json'),'closing_source_commit':commit,'closing_source_tag':tag,'records':rows,'scope':'Every currentWORKfile recovered; exact recorded births and accessible scoped consumer check pass. No hidden/global/scheduled absence claim. CI can remain running without local consumer dependency.'},indent=2)+'\n')
 print(json.dumps({'files':len(rows),'bytes':sum(r['bytes'] for r in rows),'source_commit':commit,'source_tag':tag,'ledger':seal(P/'closing-ledger.json')}))
else:
 ledger=read(P/'closing-ledger.json');tag=ledger['closing_source_tag'];remote=git('ls-remote','origin','refs/tags/'+tag)
 assert remote.split()[0]==ledger['closing_source_commit']
 files=sorted(p for p in W.rglob('*') if p.is_file());expected={r['origin']:{k:r[k] for k in ('bytes','sha256')} for r in ledger['records']}
 assert {str(p.relative_to(R)):seal(p) for p in files}==expected
 for r in ledger['records']:assert seal(R/r['recovery'])=={k:r[k] for k in ('bytes','sha256')}
 shutil.rmtree(W)
 state={'source_revision':git('rev-parse','HEAD'),'clean':not git('status','--porcelain'),'development_tracking':git('rev-list','--left-right','--count','HEAD...@{u}'),'main_tracking':git('rev-list','--left-right','--count','main...origin/main'),'worktrees':git('worktree','list'),'work_absent':not W.exists(),'out_entries':sorted(p.name for p in (R/'out').iterdir()),'deleted_files':len(files),'deleted_bytes':sum(p['bytes'] for p in expected.values()),'closing_ledger':seal(P/'closing-ledger.json'),'committed_and_pushed':True,'scope':'Finite diagnostic removal; original provider/tools/product paths preserved. Scientific and release issues remain open.'}
 assert state['clean'] and state['development_tracking']==state['main_tracking']=='0\t0'
 (P/'final-state.json').write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state))
