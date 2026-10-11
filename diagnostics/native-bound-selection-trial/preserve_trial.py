"""Preserve the rejected native trial, then remove its finished WORK and exact stage."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
TRIAL='99a4a6d3a6bff69d9f20c995383a5f2e6b8cb901';ACCEPTED='2053ef79c90752d37b0075b8137cfb9b571b6820'
ARCHIVE=ROOT/'attestations/native-vulkan/99a4a6d/native-bound-selection-trial'
NATIVE=ROOT/'bin/windows-msvc/native-99a4a6d-bound-select-trial'
TAG='evidence/issue-44-native-bound-selection-trial-99a4a6d'
PREFIX='diagnostics/native-bound-selection-trial/'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
 with p.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
 return {'bytes':p.stat().st_size,'sha256':digest}
def git(*args,env=None,input=None):return subprocess.check_output(['git',*args],env=env,input=input,text=True,timeout=30).strip()
def files(base):
 entries=list(base.rglob('*'));assert all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in entries)
 return sorted(p for p in entries if p.is_file())
def copy_all():
 records=[]
 for p in files(WORK):
  target=ARCHIVE/'files'/p.relative_to(WORK);target.parent.mkdir(parents=True,exist_ok=True)
  if target.exists():assert seal(target)==seal(p)
  else:shutil.copyfile(p,target)
  assert seal(p)==seal(target)
  records.append({'origin':str(p.relative_to(ROOT)),'recovery':str(target.relative_to(ROOT)),**seal(p)})
 return records
assert len(sys.argv)==2 and sys.argv[1] in ('archive','remove')
assert WORK==ROOT/'out/native-bound-selection-trial' and not WORK.is_symlink()
assert not git('status','--porcelain')
assert git('rev-list','--left-right','--count','HEAD...@{u}')=='0\t0'
assert git('rev-list','--left-right','--count','main...origin/main')=='0\t0'
assert git('worktree','list','--porcelain').count('worktree ')==1
if sys.argv[1]=='archive':
 assert git('rev-parse','HEAD')==TRIAL
 assert read(WORK/'retention-review.json')['accepted'] and read(WORK/'final-ownership-review.json')['accepted']
 assert read(WORK/'final-ownership-review.json')['reports']==read(WORK/'retention-review.json')['reports']
 decision=read(WORK/'root-retention-decision.json');assert len(decision['retention_checks'])==44 and len(decision['failed_conditions'])==32
 for name,expected in read(WORK/'retention-review.json')['reports'].items():assert seal(WORK/name)==expected
 review=read(WORK/'preservation-source-review.json');assert review['accepted'] and seal(Path(__file__))==review['helper_seal']
 assert not ARCHIVE.exists() and not git('tag','--list',TAG)
 useful_json={'native-cases.json','source-review.json','build-controller-review.json','product-controller-review.json','compiled-review.json','controller-review.json','analysis-review.json','controls-reader-v2-review.json','native-preparation-review.json','native-launch-review.json','retention-review.json','final-ownership-review.json','preservation-source-review.json','restore-build-controller-review.json','trial-retirement-source-review.json','restored-product-reader-review.json'}
 env=dict(os.environ,GIT_INDEX_FILE=str(WORK/'preservation.index'));blobs={}
 try:
  git('read-tree',TRIAL,env=env)
  for p in files(WORK):
   if p.suffix not in ('.py','.ps1','.patch','.slang') and p.name not in useful_json:continue
   if 'owner-temp' in p.parts:continue
   relative=str(p.relative_to(WORK));blob=git('hash-object','-w',str(p));name=PREFIX+relative
   git('update-index','--add','--cacheinfo','100644',blob,name,env=env);blobs[relative]={'git_path':name,'git_blob':blob}
  tree=git('write-tree',env=env)
  commit=git('commit-tree',tree,'-p',TRIAL,input='Preserve rejected native bound-selector trial (#44)\n\nEight original controls and eight ABBA commands completed normally.\nAll 60 complete buffers and 51 ordered markers agree; 32 of 44 frozen benefit conditions fail.\nRetain exact candidate source, reviewed controllers/readers and numerical/ownership evidence.\nReject optimization and restore accepted 2053; full rendering/qualification remains open.\n')
  git('tag',TAG,commit)
 finally:
  for p in WORK.glob('preservation.index*'):p.unlink()
 ARCHIVE.mkdir(parents=True);records=copy_all()
 for row in records:row.update(blobs.get(str(Path(row['origin']).relative_to(WORK.relative_to(ROOT))),{}))
 manifest={'source_revision':TRIAL,'restoration_revision':ACCEPTED,'source_tag':TAG,'source_commit':commit,'records':records,'scope':'Useful source and concise rejection in pushed Git tag; all raw trial files locally archived. No remote backup of raw native binary/evidence or qualification claim.'}
 (ARCHIVE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps({'source_tag':TAG,'source_commit':commit,'files':len(records),'bytes':sum(row['bytes'] for row in records)}))
else:
 assert git('rev-parse','HEAD')==ACCEPTED and git('branch','--show-current')=='development/retained-dense-step-control-2026-10-06'
 review=read(WORK/'cleanup-ready-review.json');assert review['accepted'] and review['helper_seal']==seal(Path(__file__))
 assert read(WORK/'publication-readback.json')['exact']
 assert read(WORK/'restored-source-build.json')['accepted'] and read(WORK/'restored-products.json')['accepted']
 bridge=read(WORK/'retirement-run.json');census=read(WORK/'native-retirement.json')
 assert bridge['accepted'] and not bridge['errors'] and bridge['inputs_unchanged'] and census['pass'] and not census['matches'] and not census['errors']
 assert len(census['targets'])==32 and all(t['recorded_birth_absent'] for t in census['targets'])
 for name,expected in bridge['input_seals'].items():assert seal(Path(name))==expected
 assert census['native_stage']==subprocess.check_output(['wslpath','-w',str(NATIVE)],text=True,timeout=10).strip()
 assert census['stage']==subprocess.check_output(['wslpath','-w',str(WORK)],text=True,timeout=10).strip()
 targets=[]
 for action in read(WORK/'native-cases.json')['actions']:
  owner=read(WORK/action/'run/owner.json')
  targets.extend([(action+' owner',int(owner['owner_pid']),owner['owner_start_time_utc']),(action+' child',int(owner['test_pid']),owner['test_start_time_utc'])])
 assert len(set(targets))==32 and [(t['kind'],int(t['pid']),t['birth']) for t in census['targets']]==targets
 manifest=read(ARCHIVE/'manifest.json');assert git('ls-remote','origin','refs/tags/'+TAG).split()[0]==manifest['source_commit']
 for row in manifest['records']:assert seal(ROOT/row['origin'])==seal(ROOT/row['recovery'])=={k:row[k] for k in ('bytes','sha256')}
 plan=read(WORK/'execution-plan.json');prefix=str(NATIVE.relative_to(ROOT))+'/'
 expected={name:record for name,record in plan['input_seals'].items() if name.startswith(prefix)}
 assert {str(p.relative_to(ROOT)):seal(p) for p in files(NATIVE)}==expected
 retained={(row['bytes'],row['sha256']) for row in manifest['records']}
 assert all((row['bytes'],row['sha256']) in retained for row in expected.values()),'Native stage file has no verified retained archive counterpart'
 births=[]
 def absent(pid,ticks):
  try:actual=Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
  except (FileNotFoundError,ProcessLookupError):actual=None
  assert ticks and actual!=ticks
  births.append({'pid':pid,'start_ticks':ticks,'recorded_birth_absent':True})
 for action in read(WORK/'native-cases.json')['actions']:
  b=read(WORK/action/'native-control-bridge.json');rr=read(WORK/action/'retirement-run.json')
  absent(b['linux_owner_pid'],b['linux_owner_start_ticks']);absent(b['linux_child_pid'],b['linux_child_start_ticks']);absent(rr['pid'],rr['start_ticks'])
 matches=[];hidden=[];visible=0
 for p in Path('/proc').iterdir():
  if not p.name.isdigit() or int(p.name) in (os.getpid(),os.getppid()):continue
  try:command=(p/'cmdline').read_bytes().replace(b'\0',b' ');visible+=1
  except (FileNotFoundError,ProcessLookupError,PermissionError) as error:hidden.append({'pid':int(p.name),'reason':type(error).__name__});continue
  if any(value.encode() in command for value in (str(WORK),str(WORK.relative_to(ROOT)),str(NATIVE),str(NATIVE.relative_to(ROOT)))):matches.append({'pid':int(p.name),'command':command.decode(errors='replace')})
 assert not matches,matches
 records=copy_all();publication=ARCHIVE/'publication';publication.mkdir()
 ledger={'source_tag':TAG,'source_commit':manifest['source_commit'],'manifest_seal':seal(ARCHIVE/'manifest.json'),'records':records,'native_stage_files':len(expected),'native_stage_bytes':sum(v['bytes'] for v in expected.values()),'linux_births':births,'linux_visible_processes':visible,'linux_hidden_records':hidden,'scope':'Fresh accessible consumer scan and exact saved births; no global/hidden/scheduled absence claim.'}
 (publication/'closing-ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
 shutil.rmtree(NATIVE);shutil.rmtree(WORK)
 state={'source_revision':ACCEPTED,'clean':not git('status','--porcelain'),'tracking':git('rev-list','--left-right','--count','HEAD...@{u}'),'main_tracking':git('rev-list','--left-right','--count','main...origin/main'),'worktrees':git('worktree','list'),'finished_work_absent':not WORK.exists(),'finished_native_stage_absent':not NATIVE.exists(),'out_entries':sorted(p.name for p in (ROOT/'out').iterdir()),'removed_work_files':len(records),'removed_work_bytes':sum(v['bytes'] for v in records),'scope':'Rejected trial source/raw evidence preserved; acceptedsource/products restored. Fullscientific/rendering/release qualification remains open.'}
 assert state['clean'] and state['tracking']==state['main_tracking']=='0\t0'
 (publication/'removal-state.json').write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state))
