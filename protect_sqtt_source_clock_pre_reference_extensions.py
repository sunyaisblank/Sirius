"""Push useful source, preserve exact consumed versions, then close only this batch."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,zipfile
ROOT=Path.cwd();W=Path(__file__).resolve().parent
assert ROOT==Path('/home/astra/.project/Sirius') and W==ROOT/'out/native-sqtt-source-clock'
assert len(sys.argv)==2 and sys.argv[1] in ('source','archive','review','remove')
mode=sys.argv[1]
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def expected(v):return {k:v[k] for k in ('bytes','sha256')}
def key(v):return v['bytes'],v['sha256']
def git(*a,**kwargs):return subprocess.check_output(['git',*a],**kwargs).decode().strip()
head=read(W/'plan.json')['source_revision'];D=ROOT/'attestations/native-vulkan'/head[:7]/'native-sqtt-source-clock';tag='evidence/issue-44-native-sqtt-source-clock-'+head[:7]
def source_check():
 assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}')
 assert not git('status','--porcelain')
 assert git('rev-parse','HEAD^{tree}')==git('rev-parse','b81e061^{tree}')
 assert git('rev-parse','main')==git('rev-parse','origin/main')=='5e706073a3fd2f9afcf8efbcde2f094122050f6a'
 assert len([v for v in git('worktree','list','--porcelain').splitlines() if v.startswith('worktree ')])==1
source_check()
review=read(W/'preservation-source-review.json');assert review['accepted'] and seal(Path(__file__))==review['controller_seal']
assert read(W/'batch-disposition.json')['finished_offline_batch']
if mode=='source':
 assert not git('ls-remote','origin','refs/tags/'+tag)
 assert subprocess.run(['git','show-ref','--verify','--quiet','refs/tags/'+tag]).returncode!=0
 paths=sorted(p for p in W.iterdir() if p.is_file() and (p.suffix in ('.py','.cpp','.cc','.hpp','.h','.md','.txt','.json') or p.name.endswith('-LICENSE')) and p.name not in ('saved-isa-map.json','saved-isa2-map.json','isa-result.json','isa2-result.json','issue-44-before.json'))
 assert any(p.name=='decoder-LICENSE' for p in paths) and W/'saved-isa-summary.json' in paths
 index=W/'preservation.index';assert not index.exists() and not (W/'source-preservation.json').exists()
 env=dict(os.environ,GIT_INDEX_FILE=str(index));subprocess.run(['git','read-tree','--empty'],env=env,check=True);records=[]
 for p in paths:
  assert not p.is_symlink();data=p.read_bytes();blob=git('hash-object','-w','--stdin',input=data)
  subprocess.run(['git','update-index','--add','--cacheinfo','100644,'+blob+','+p.name],env=env,check=True)
  records.append({'origin':str(p.relative_to(ROOT)),'path':p.name,'blob':blob,**seal(p)})
 tree=git('write-tree',env=env);commit=git('commit-tree',tree,input=('Preserve reviewed saved SQTT source and event-clock joins at '+head[:7]+' (#44)\n\nNo product change or scientific/release qualification.\n').encode())
 subprocess.run(['git','tag',tag,commit],check=True);subprocess.run(['git','push','origin','refs/tags/'+tag],check=True)
 assert git('ls-remote','origin','refs/tags/'+tag).split()[0]==commit;index.unlink();source_check()
 (W/'source-preservation.json').write_text(json.dumps({'tag':tag,'commit':commit,'source_revision':head,'pushed':True,'remote_commit_verified':True,'files':records},indent=2)+'\n')
 print(json.dumps({'tag':tag,'commit':commit,'files':len(records),'pushed':True}));sys.exit(0)
source=read(W/'source-preservation.json');assert source['pushed'] and source['tag']==tag and git('rev-parse',tag)==source['commit']
assert git('ls-remote','origin','refs/tags/'+tag).split()[0]==source['commit']
recoveries={};verified=set();zips={}
def verify(value,recovery):
 marker=key(value),json.dumps(recovery,sort_keys=True)
 if marker in verified:return
 if recovery['kind']=='immutable-Git-blob':
  assert git('rev-parse',recovery['revision']+':'+recovery['path'])==recovery['blob']
  data=subprocess.check_output(['git','cat-file','blob',recovery['blob']]);assert {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}==expected(value)
 elif recovery['kind']=='protected-ZIP-member':
  p=ROOT/recovery['archive_path']
  if str(p) not in zips:assert seal(p)==recovery['archive_seal'];zips[str(p)]=zipfile.ZipFile(p)
  z=zips[str(p)];info=z.getinfo(recovery['member']);assert info.file_size==value['bytes'] and info.CRC==recovery['CRC32']
  data=z.read(info);assert {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}==expected(value)
 else:
  assert recovery['kind'] in ('local-archive-payload','local-publication-payload','retained-prior-payload','retained-native-export','retained-provider-payload','retained-archive-payload')
  p=ROOT/recovery['path'];assert p.is_file() and not p.is_symlink() and seal(p)==expected(value)
 verified.add(marker);recoveries[key(value)]=recovery
for v in source['files']:recoveries[key(v)]={'kind':'immutable-Git-blob','revision':source['commit'],'blob':v['blob'],'path':v['path']}
if mode=='archive':
 assert not D.exists()
 priors=[]
 for name in ('native-spm-cost','native-counter-views','native-counter-clock'):
  p=ROOT/'attestations/native-vulkan'/head[:7]/name/'manifest.json';m=read(p);priors.append({'path':str(p.relative_to(ROOT)),**seal(p)})
  for v in [*m['executed_input_ledger'],*m['records'],*([m['tool_archive']] if 'tool_archive' in m else [])]:recoveries[key(v)]=v['recovery']
 # Current snapshot takes precedence over ancestor source recovery aliases.
 for v in source['files']:recoveries[key(v)]={'kind':'immutable-Git-blob','revision':source['commit'],'blob':v['blob'],'path':v['path']}
 payloads={}
 def resolve(p,v):
  if key(v) not in recoveries:
   assert p.is_file() and not p.is_symlink() and seal(p)==expected(v),str(p)
   target=D/'payloads'/v['sha256'];target.parent.mkdir(parents=True,exist_ok=True);assert not target.exists();shutil.copy2(p,target);assert seal(target)==expected(v)
   payloads[str(target.relative_to(D))]=expected(v);recoveries[key(v)]={'kind':'local-archive-payload','path':str(target.relative_to(ROOT))}
  r=recoveries[key(v)];verify(v,r);return r
 ledger={}
 def collect(origin,v,provenance):
  identity=origin,*key(v);p=Path(origin) if Path(origin).is_absolute() else ROOT/origin
  if identity not in ledger:ledger[identity]={'origin':origin,**expected(v),'ledgers':[],'recovery':resolve(p,v)}
  ledger[identity]['ledgers'].append(provenance)
 for p in sorted(W.glob('*.json')):
  def visit(node):
   if isinstance(node,dict):
    for name,child in node.items():
     if name in ('input_seals','source_input_seals','files') and isinstance(child,dict):
      for origin,v in child.items():
       if isinstance(v,dict) and {'bytes','sha256'}<=v.keys():collect(origin,v,str(p.relative_to(ROOT)))
     else:visit(child)
   elif isinstance(node,list):
    for child in node:visit(child)
  visit(read(p))
 tool=ROOT/read(W/'cpu-decoder-preparation.json')['extracted_library']['path'];records=[]
 for p in [*sorted(W.rglob('*')),tool]:
  assert not p.is_symlink()
  if p.is_file():v=seal(p);records.append({'origin':str(p.relative_to(ROOT)),**v,'recovery':resolve(p,v)})
 manifest={'source_revision':head,'source_snapshot':source,'disposition':read(W/'batch-disposition.json'),'prior_manifests':priors,'executed_input_ledger':list(ledger.values()),'records':records,'payloads':payloads,'scope':'Exact finite offline SQTT code/source and event-clock observations. Useful source pushed; official binary package/library and detailed raw observations local. No product change, optimization, qualification or cleanup verdict.'}
 D.mkdir(parents=True,exist_ok=True);(D/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 receipt={'archive':str(D.relative_to(ROOT)),'manifest':seal(D/'manifest.json'),'records':len(records),'executed_origin_version_pairs':len(ledger),'payload_count':len(payloads),'payload_bytes':sum(v['bytes'] for v in payloads.values()),'all_recovery_bytes_verified':True}
 (W/'protection.json').write_text(json.dumps(receipt,indent=2)+'\n');source_check();print(json.dumps(receipt));sys.exit(0)
protection=read(W/'protection.json');assert protection['archive']==str(D.relative_to(ROOT)) and seal(D/'manifest.json')==protection['manifest'];manifest=read(D/'manifest.json');assert manifest['source_snapshot']==source
for v in [*manifest['executed_input_ledger'],*manifest['records']]:verify(v,v['recovery'])
for name,v in manifest['payloads'].items():assert seal(D/name)==v
for v in source['files']:assert seal(ROOT/v['origin'])==expected(v)
for z in zips.values():z.close()
if mode=='review':
 p=W/'preservation-reviewed.json';assert not p.exists();p.write_text(json.dumps({'accepted':True,'manifest':protection['manifest'],'source_tag':tag,'recovery_versions':len(verified),'scope':'Root executed whole-byte recovery verification only; independent reviewer must accept actual manifest before removal.'},indent=2)+'\n');print(json.dumps({'accepted':True,'recovery_versions':len(verified)}));sys.exit(0)
assert read(W/'preservation-reviewed.json')['manifest']==protection['manifest']
independent=read(W/'preservation-manifest-review.json');assert independent['accepted'] and independent['manifest']==protection['manifest']
for name in ('dump-run.json','isa-run.json','isa2-run.json'):
 run=read(W/name);assert run['returncode']==(1 if name=='dump-run.json' else 0) and not run['errors'] and not run['deferred_signals'] and not run['limit_reason'] and run['final_owned_group']==[]
assert 'undefined symbol: rocprof_trace_decoder_dump_data' in (W/'dump-stderr.txt').read_text()
assert not (W/'dump-result.json').exists() and not list(W.glob('dump-tokens-*.bin'))
def birth(pid):
 try:
  s=Path('/proc',str(pid),'stat').read_text();f=s[s.rfind(')')+2:].split();return int(f[19])
 except (FileNotFoundError,ProcessLookupError):return None
known={}
for p in W.glob('*.json'):
 def visit_births(v):
  if isinstance(v,dict):
   if isinstance(v.get('pid'),int) and isinstance(v.get('start_ticks'),(str,int)):known[(v['pid'],int(v['start_ticks']))]={'pid':v['pid'],'start_ticks':int(v['start_ticks'])}
   for c in v.values():visit_births(c)
  elif isinstance(v,list):
   for c in v:visit_births(c)
 visit_births(read(p))
for v in known.values():assert birth(v['pid'])!=v['start_ticks'],v
def tree(folder):
 assert folder.is_dir() and not folder.is_symlink();result={}
 for p in sorted(folder.rglob('*')):
  assert not p.is_symlink() and (p.is_file() or p.is_dir());result[str(p.relative_to(ROOT))]=seal(p) if p.is_file() else 'directory'
 return result
retain=['bin/linux-gcc/qualification-tools/rocprof-trace-decoder-0.1.3','bin/windows-msvc/native-b81e061-executed-cost','bin/windows-msvc/qualification-tools/rdts-2026-10-02-1872','bin/linux-gcc/src/sirius/backend/retained','bin/linux-gcc/tests/backend/resources','bin/linux-gcc/src/sirius/app/resources']
retained={name:tree(ROOT/name) for name in retain};zip_path=ROOT/'attestations/native-vulkan/6be3132/native-executed-cost/payloads/128fdebb889d1fda1d04ba9e4792fbef1b18c086b4e901642262653a7a53917e';zip_seal=seal(zip_path);assert zip_seal=={'bytes':424481206,'sha256':'128fdebb889d1fda1d04ba9e4792fbef1b18c086b4e901642262653a7a53917e'}
assert [sum(isinstance(v,dict) for v in t.values()) for t in retained.values()]==[1,94,70,35,23,32]
closing=D/'closing';assert not closing.exists();before=tree(W);entries=[]
for origin,v in before.items():
 if not isinstance(v,dict):continue
 r=recoveries.get(key(v))
 if r is None:
  p=closing/'payloads'/v['sha256'];p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/origin,p);r={'kind':'local-publication-payload','path':str(p.relative_to(ROOT))}
 verify(v,r);entries.append({'origin':origin,**v,'recovery':r})
matches=[];errors=[];accessible=0
needles=('out/native-sqtt-source-clock','native-b81e061-executed-cost','rdts-2026-10-02-1872','qualification-tools/rocprof-trace-decoder-0.1.3')
for p in Path('/proc').glob('[0-9]*'):
 if int(p.name)==os.getpid():continue
 accessible+=1
 for field in ('cmdline','environ','maps','cwd','exe'):
  try:
   value=os.readlink(p/field) if field in ('cwd','exe') else (p/field).read_bytes().decode(errors='replace')
   if any(n in value for n in needles):matches.append({'pid':int(p.name),'field':field})
  except (OSError,ValueError):errors.append({'pid':int(p.name),'field':field})
assert not matches,matches
assert tree(W)==before;source_check()
closing.mkdir(parents=True,exist_ok=True);(closing/'manifest.json').write_text(json.dumps({'source_revision':head,'source_snapshot':source,'deletion_recovery_ledger':entries,'retained_trees':retained,'protected_zip':{'path':str(zip_path.relative_to(ROOT)),**zip_seal}},indent=2)+'\n')
receipt={'source_revision':head,'source_committed_pushed':True,'source_tag':tag,'manifest':protection['manifest'],'closing_manifest':seal(closing/'manifest.json'),'remove':['out/native-sqtt-source-clock'],'retained':retain,'recorded_births':list(known.values()),'recorded_births_absent':True,'closing_owner':{'pid':os.getpid(),'start_ticks':birth(os.getpid())},'accessible_linux_processes':accessible,'consumer_matches':matches,'consumer_read_errors':errors,'removed':False,'files':len(entries),'bytes':sum(v['bytes'] for v in entries),'scope':'Scoped saved Linux births and accessible consumer fields only. No hidden/global/scheduled absence, Windows census or scientific/release qualification claim.'}
(D/'cleanup.json').write_text(json.dumps(receipt,indent=2)+'\n')
assert tree(W)==before and seal(zip_path)==zip_seal
for name,v in retained.items():assert tree(ROOT/name)==v
for v in known.values():assert birth(v['pid'])!=v['start_ticks']
source_check();shutil.rmtree(W);assert not W.exists() and seal(zip_path)==zip_seal
for name,v in retained.items():assert tree(ROOT/name)==v
source_check();receipt.update(removed=True,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),retained_verified_after=True,out_empty=not any((ROOT/'out').iterdir()),final_tracking=git('rev-list','--left-right','--count','HEAD...@{u}'))
(D/'cleanup.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ('removed','files','bytes','out_empty','final_tracking')}))
