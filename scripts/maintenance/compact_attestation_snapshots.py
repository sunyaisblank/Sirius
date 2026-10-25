#!/usr/bin/env python3
"""Share identical immutable historical snapshots without changing their paths.

Only untracked, inactive attestation files with equal content and metadata are
eligible. Source/build trees, archives, forensic originals, and the current
qualification/dependency closure are excluded. No payload archive is created.
Use restore to give a linked historical snapshot an independent inode before
ever editing it. Nothing here schedules automatic cleanup or rewrites Git.
"""
from pathlib import Path
from collections import defaultdict
import argparse,base64,hashlib,json,os,shutil,stat,subprocess,time

R=Path(__file__).resolve().parents[2]
ROOT=R/'attestations'
SCRIPT=Path(__file__).resolve()
EXCLUDED={'qualification','forensic-images','forensic','cleanup'}
ARCHIVES=('.zip','.gz','.bz2','.xz','.zst','.7z','.tar','.tgz','.iso','.vhd','.vhdx')

def digest(path):
 with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def metadata(path):
 s=path.lstat();assert stat.S_ISREG(s.st_mode) and not path.is_symlink()
 attrs={n:base64.b64encode(os.getxattr(path,n,follow_symlinks=False)).decode() for n in sorted(os.listxattr(path,follow_symlinks=False))}
 return dict(bytes=s.st_size,mode=stat.S_IMODE(s.st_mode),mtime_ns=s.st_mtime_ns,uid=s.st_uid,gid=s.st_gid,device=s.st_dev,inode=s.st_ino,links=s.st_nlink,allocated_bytes=s.st_blocks*512,xattrs_sha256=hashlib.sha256(json.dumps(attrs,sort_keys=True).encode()).hexdigest())

def content_metadata(value):
 return tuple(value[k] for k in ('bytes','mode','mtime_ns','uid','gid','device','xattrs_sha256'))

def stat_identity(value):return tuple(value[k] for k in ('bytes','mode','mtime_ns','uid','gid','device','inode','xattrs_sha256'))

def relative(path):return str(path.relative_to(R))

def target(name):
 p=R/name
 assert p.is_relative_to(ROOT) and not any(part=='.git' for part in p.parts)
 assert p.resolve(strict=True)==p and p.is_file() and not p.is_symlink()
 return p

def write(path,value):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')

def query(*args):return subprocess.check_output(args,cwd=R,text=True)

def protected():
 paths=set()
 for root in ('native-cpu-local-identities','native-cpu-local-identities-cohort','native-cpu-local-identities-compatibility','native-cpu-owned-lower-increments','native-cpu-coherent-profile'):
  for p in (R/'out'/root).rglob('actions.json'):
   for action in json.loads(p.read_text()).values():paths.update(str((R/n).resolve()) for n in action.get('frozen_inputs',[]))
 return paths

def busy():
 """Read only same-user process metadata; return only paths in this checkout."""
 found=set();checked=0
 code='''from pathlib import Path
import json,os,sys
p=Path(sys.argv[1]);root=sys.argv[2];found=set()
try:
 if p.stat().st_uid==int(sys.argv[3]):
  for line in (p/'maps').read_text().splitlines():
   fields=line.split(None,5)
   if len(fields)==6 and fields[5].startswith(root+'/'):found.add(fields[5])
  for fd in (p/'fd').iterdir():
   try:q=os.readlink(fd)
   except FileNotFoundError:continue
   if q.startswith(root+'/'):found.add(q)
except (FileNotFoundError,ProcessLookupError):pass
print(json.dumps(sorted(found)))
'''
 for p in Path('/proc').iterdir():
  if not p.name.isdecimal():continue
  try:
   if p.stat().st_uid!=os.getuid():continue
   try:
    for line in (p/'maps').read_text().splitlines():
     fields=line.split(None,5)
     if len(fields)==6 and fields[5].startswith(str(R)+'/'):found.add(fields[5])
    for fd in (p/'fd').iterdir():
     try:q=os.readlink(fd)
     except FileNotFoundError:continue
     if q.startswith(str(R)+'/'):found.add(q)
   except PermissionError:
    result=subprocess.check_output(['sudo','-n','/usr/bin/python3','-B','-c',code,str(p),str(R),str(os.getuid())],text=True)
    found.update(json.loads(result))
   checked+=1
  except (FileNotFoundError,ProcessLookupError):continue
 return found,checked

def make_plan(destination):
 assert not destination.exists()
 tracked=set(subprocess.check_output(['git','ls-files','-z'],cwd=R).split(b'\0'))
 excluded=protected();opened,processes=busy();excluded|=opened
 groups=defaultdict(list);scanned=0
 for base,dirs,files in os.walk(ROOT,followlinks=False):
  root=Path(base);parts=root.relative_to(ROOT).parts
  dirs[:]=[d for d in dirs if d!='.git' and not (root/d).is_symlink()]
  if parts and parts[0] in EXCLUDED:dirs[:]=[];continue
  for name in files:
   p=root/name
   if p.is_symlink() or str(p) in excluded or os.fsencode(relative(p)) in tracked or name.lower().endswith(ARCHIVES):continue
   s=p.lstat()
   if not stat.S_ISREG(s.st_mode) or s.st_size<1048576 or s.st_uid!=os.getuid():continue
   m=metadata(p);groups[content_metadata(m)].append((p,m));scanned+=1
 actions=[];hashed=0;integrity=[]
 for _,rows in sorted(groups.items()):
  if len({m['inode'] for _,m in rows})<2:continue
  content=defaultdict(list)
  for p,before in rows:
   sha=digest(p);after=metadata(p);assert stat_identity(before)==stat_identity(after),p
   if len(p.name)==64 and all(c in '0123456789abcdef' for c in p.name) and p.name!=sha:
    integrity.append(relative(p));continue
   content[sha].append((p,before));hashed+=1
   if hashed%200==0:print(json.dumps(dict(hashed_files=hashed,planned_aliases=len(actions))),flush=True)
  for sha,items in content.items():
   if len({m['inode'] for _,m in items})<2:continue
   items.sort(key=lambda item:(item[1]['allocated_bytes'],relative(item[0])))
   canonical,cm=items[0]
   seen={cm['inode']}
   for p,m in items[1:]:
    if m['inode'] in seen:continue
    seen.add(m['inode'])
    if m['links']!=1:continue
    actions.append(dict(canonical=relative(canonical),path=relative(p),sha256=sha,canonical_before=cm,before=m))
 plan=dict(schema=1,root=str(R),created_unix=time.time(),source=query('git','rev-parse','HEAD').strip(),script_sha256=digest(SCRIPT),scanned_large_files=scanned,hashed_files=hashed,same_user_processes_checked=processes,protected_paths=len(excluded),actions=actions,gross_allocated_bytes_reclaimable=sum(a['before']['allocated_bytes'] for a in actions),checksum_name_mismatches=integrity,
  scope='Physical sharing of immutable historical attestation snapshots only. Every payload independently SHA256 verified; size/mode/mtime/owner/group/device/xattrs equal. Existing paths/content/required metadata remain. Current frozen inputs, mapped/open files, tracked files, archives and original forensic/qualification evidence excluded. No source/build-tree move, deletion of unique data, Git history rewrite or new binary archive. Restore makes a path independent again before any edit.')
 write(destination,plan);print(json.dumps({k:plan[k] for k in ('scanned_large_files','hashed_files','gross_allocated_bytes_reclaimable','checksum_name_mismatches')}|dict(planned_aliases=len(actions))))

def apply_plan(path,limit):
 plan=json.loads(path.read_text());assert plan['root']==str(R) and plan['script_sha256']==digest(SCRIPT)
 journal=path.with_name('operations.jsonl');opened,processes=busy();excluded=protected()|opened;done=0;reclaimed=0;already=0
 before_free=shutil.disk_usage(R).free
 with journal.open('a',buffering=1) as log:
  for row in plan['actions']:
   a,b=target(row['canonical']),target(row['path']);am,bm=metadata(a),metadata(b)
   assert content_metadata(am)==content_metadata(row['canonical_before']) and content_metadata(bm)==content_metadata(row['before'])
   if a.samefile(b):already+=1;continue
   if limit and done>=limit:break
   assert str(a) not in excluded and str(b) not in excluded
   assert stat_identity(am)==stat_identity(row['canonical_before']) and stat_identity(bm)==stat_identity(row['before']) and bm['links']==1
   assert digest(a)==digest(b)==row['sha256']
   assert stat_identity(metadata(a))==stat_identity(am) and stat_identity(metadata(b))==stat_identity(bm)
   temp=b.with_name('.sirius-dedup-'+str(os.getpid()));assert not temp.exists()
   log.write(json.dumps(dict(status='prepared',path=row['path'],canonical=row['canonical'],sha256=row['sha256']))+'\n');log.flush();os.fsync(log.fileno())
   try:
    os.link(a,temp,follow_symlinks=False)
    assert temp.samefile(a) and stat_identity(metadata(b))==stat_identity(bm)
    os.replace(temp,b)
   finally:
    if temp.exists():temp.unlink()
   assert a.samefile(b) and content_metadata(metadata(b))==content_metadata(row['before']) and digest(b)==row['sha256']
   done+=1;reclaimed+=bm['allocated_bytes'];log.write(json.dumps(dict(status='complete',path=row['path'],canonical=row['canonical'],sha256=row['sha256'],gross_allocated_bytes_retired=bm['allocated_bytes']))+'\n')
   if done%100==0:print(json.dumps(dict(new_aliases=done,gross_allocated_bytes_retired=reclaimed)),flush=True)
  log.flush();os.fsync(log.fileno())
 summary=dict(time_unix=time.time(),source=query('git','rev-parse','HEAD').strip(),new_aliases=done,already_shared=already,gross_allocated_bytes_retired=reclaimed,filesystem_free_bytes_before=before_free,filesystem_free_bytes_after=shutil.disk_usage(R).free,same_user_processes_checked=processes)
 name='batch-'+str(time.time_ns())+'.json';write(path.parent/name,summary);print(json.dumps(summary))

def verify(path):
 plan=json.loads(path.read_text());done=0;remaining=0
 for row in plan['actions']:
  a,b=target(row['canonical']),target(row['path']);am,bm=metadata(a),metadata(b)
  assert content_metadata(am)==content_metadata(row['canonical_before']) and content_metadata(bm)==content_metadata(row['before'])
  if a.samefile(b):assert digest(b)==row['sha256'];done+=1
  else:assert stat_identity(bm)==stat_identity(row['before']);remaining+=1
 result=dict(verified_shared_paths=done,remaining_independent_paths=remaining,plan_sha256=digest(path),all_paths_content_and_metadata_valid=True)
 print(json.dumps(result));return result

def restore(path,name):
 plan=json.loads(path.read_text());rows=[r for r in plan['actions'] if r['path']==name];assert len(rows)==1;row=rows[0]
 p=target(name);assert digest(p)==row['sha256'];before=metadata(p)
 temp=p.with_name('.sirius-independent-'+str(os.getpid()));assert not temp.exists()
 try:
  shutil.copy2(p,temp,follow_symlinks=False);assert digest(temp)==row['sha256'] and content_metadata(metadata(temp))==content_metadata(before)
  os.replace(temp,p)
 finally:
  if temp.exists():temp.unlink()
 assert metadata(p)['links']==1 and digest(p)==row['sha256'];print(json.dumps(dict(restored_independent_path=name)))

parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['plan','apply','verify','restore']);parser.add_argument('plan',type=Path);parser.add_argument('--limit',type=int,default=0);parser.add_argument('--path')
args=parser.parse_args();plan=args.plan.resolve();assert plan.is_relative_to(R/'attestations/cleanup')
if args.mode=='plan':make_plan(plan)
elif args.mode=='apply':apply_plan(plan,args.limit)
elif args.mode=='verify':verify(plan)
else:assert args.path;restore(plan,args.path)
