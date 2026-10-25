#!/usr/bin/env python3
"""Archive inactive historical payloads once per SHA256, then retire expansions.

Metadata, source, current inputs, Git-tracked files, forensic originals and
existing recovery archives remain in place. `restore PLAN --path OLD_PATH`
restores an exact payload without overwriting an existing file. Compressed
objects remain local under the plan directory; the index is suitable for Git.
"""
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,gzip,hashlib,json,os,shutil,stat,subprocess,time

R=Path(__file__).resolve().parents[2]
BASE=R/'scripts/maintenance/compact_attestation_snapshots.py'
namespace={'__file__':str(BASE)}
exec(compile(BASE.read_text().split('parser=argparse.ArgumentParser',1)[0],str(BASE),'exec'),namespace)
metadata=namespace['metadata'];identity=namespace['stat_identity'];meta=namespace['content_metadata'];digest=namespace['digest'];busy=namespace['busy'];protected=namespace['protected'];target=namespace['target']
SCRIPT=Path(__file__).resolve()
SOURCE_SUFFIXES={'.json','.jsonl','.xml','.txt','.md','.py','.cpp','.c','.cc','.h','.hpp','.cs','.ispc','.isph','.slang','.spvasm','.sh','.ps1','.cmake','.yml','.yaml'}
ARCHIVES=namespace['ARCHIVES']
EXCLUDED={'qualification','forensic-images','forensic','cleanup'}
EMPTY_ATTRS=hashlib.sha256(b'{}').hexdigest()

def write(path,value):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')

def load(path):return json.loads(path.read_text())

def check_plan(path):
 p=load(path);assert p['root']==str(R) and p['script_sha256']==digest(SCRIPT) and p['helper_sha256']==digest(BASE)
 return p

def make_plan(path):
 assert not path.exists();tracked=set(subprocess.check_output(['git','ls-files','-z'],cwd=R).split(b'\0'))
 opened,processes=busy();excluded=protected()|opened;candidates=[]
 for base,dirs,files in os.walk(R/'attestations',followlinks=False):
  root=Path(base);parts=root.relative_to(R/'attestations').parts
  dirs[:]=[d for d in dirs if d!='.git' and not (root/d).is_symlink()]
  if parts and parts[0] in EXCLUDED:dirs[:]=[];continue
  for name in files:
   p=root/name
   if p.is_symlink() or str(p) in excluded or os.fsencode(str(p.relative_to(R))) in tracked or p.suffix.lower() in SOURCE_SUFFIXES or name.lower().endswith(ARCHIVES):continue
   s=p.lstat()
   if not stat.S_ISREG(s.st_mode) or s.st_size<1048576 or s.st_uid!=os.getuid():continue
   m=metadata(p)
   if m['xattrs_sha256']!=EMPTY_ATTRS:continue
   candidates.append((p,m))
 def inspect(item):
  p,before=item;sha=digest(p);assert identity(metadata(p))==identity(before)
  if len(p.name)==64 and all(c in '0123456789abcdef' for c in p.name):assert sha==p.name,p
  return dict(path=str(p.relative_to(R)),sha256=sha,before=before)
 records=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(inspect,item) for item in candidates]
  for future in as_completed(futures):
   records.append(future.result())
   if len(records)%500==0:print(json.dumps(dict(hashed_paths=len(records),total=len(candidates))),flush=True)
 records.sort(key=lambda r:r['path']);objects={};inodes={}
 for row in records:
  objects.setdefault(row['sha256'],dict(source_path=row['path'],bytes=row['before']['bytes']))
  inodes[(row['before']['device'],row['before']['inode'])]=row['before']
 p=dict(schema=1,root=str(R),source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),script_sha256=digest(SCRIPT),helper_sha256=digest(BASE),created_unix=time.time(),same_user_processes_checked=processes,protected_paths=len(excluded),records=records,objects=objects,expanded_unique_allocated_bytes=sum(v['allocated_bytes'] for v in inodes.values()),
  scope='Cold untracked historical binary/input/resource snapshots only, at least1MiB. Current frozen inputs and open/mapped files, Git-tracked files, source/metadata, current qualification, forensic originals and existing archives excluded. Preserve every unique byte and original file mode/mtime/owner/group in a local compressed content-addressed archive; remove only verified expanded copies. All old paths are recoverable with this script; no overwrite of newer state, Git history rewrite or system capability change.',restore_command='python3 -B scripts/maintenance/archive_cold_payloads.py restore '+str(path.relative_to(R))+' --path <original-relative-path>')
 write(path,p);print(json.dumps(dict(planned_paths=len(records),unique_payloads=len(objects),expanded_unique_allocated_bytes=p['expanded_unique_allocated_bytes'])))

def pack(path):
 p=check_plan(path);out=path.parent/'objects.json';assert not out.exists();objects={}
 def one(item):
  sha,entry=item;source=target(entry['source_path']);m=metadata(source)
  expected=next(r['before'] for r in p['records'] if r['path']==entry['source_path']);assert identity(m)==identity(expected)
  dest=path.parent/'payloads'/sha[:2]/(sha+'.gz');dest.parent.mkdir(parents=True,exist_ok=True)
  if not dest.exists():
   temp=dest.with_suffix('.partial');assert not temp.exists();observed=hashlib.sha256();size=0
   with source.open('rb') as incoming,temp.open('xb') as stream:
    with gzip.GzipFile(filename='',mode='wb',fileobj=stream,compresslevel=1,mtime=0) as compressed:
     while data:=incoming.read(1048576):observed.update(data);size+=len(data);compressed.write(data)
    stream.flush();os.fsync(stream.fileno())
   assert observed.hexdigest()==sha and size==entry['bytes'] and identity(metadata(source))==identity(m)
   os.replace(temp,dest)
  observed=hashlib.sha256();size=0
  with gzip.open(dest,'rb') as stream:
   while data:=stream.read(1048576):observed.update(data);size+=len(data)
  assert observed.hexdigest()==sha and size==entry['bytes']
  return sha,dict(path=str(dest.relative_to(R)),bytes=dest.stat().st_size,allocated_bytes=dest.stat().st_blocks*512,sha256=digest(dest),uncompressed_bytes=size,uncompressed_sha256=sha,verified=True)
 with ThreadPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(one,item) for item in p['objects'].items()]
  for f in as_completed(futures):
   sha,row=f.result();objects[sha]=row
   if len(objects)%100==0:print(json.dumps(dict(verified_objects=len(objects),total=len(p['objects']))),flush=True)
 write(out,dict(plan_sha256=digest(path),objects=objects,all_payloads_roundtrip_verified=True,compressed_allocated_bytes=sum(o['allocated_bytes'] for o in objects.values())))
 print(json.dumps(dict(verified_objects=len(objects),compressed_allocated_bytes=sum(o['allocated_bytes'] for o in objects.values()),expanded_unique_allocated_bytes=p['expanded_unique_allocated_bytes'])))

def retire(path):
 p=check_plan(path);objects=load(path.parent/'objects.json');assert objects['plan_sha256']==digest(path) and objects['all_payloads_roundtrip_verified']
 for row in objects['objects'].values():
  blob=R/row['path'];assert blob.is_file() and not blob.is_symlink() and digest(blob)==row['sha256']
 opened,processes=busy();excluded=protected()|opened
 assert all(str(R/r['path']) not in excluded for r in p['records'])
 journal=path.parent/'retired.jsonl';free_before=shutil.disk_usage(R).free;removed=0;freed=0
 with journal.open('a',buffering=1) as log:
  for row in p['records']:
   file=R/row['path']
   if not file.exists():continue
   file=target(row['path']);m=metadata(file);assert identity(m)==identity(row['before']) and digest(file)==row['sha256']
   assert identity(metadata(file))==identity(m)
   if m['links']==1:freed+=m['allocated_bytes']
   file.unlink();removed+=1;log.write(json.dumps(dict(path=row['path'],sha256=row['sha256'],archive=objects['objects'][row['sha256']]['path']))+'\n')
   if removed%100==0:log.flush();os.fsync(log.fileno())
   if removed%500==0:print(json.dumps(dict(retired_paths=removed,gross_unique_allocated_bytes_retired=freed)),flush=True)
  log.flush();os.fsync(log.fileno())
 summary=dict(source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),retired_paths_this_run=removed,total_paths_absent=sum(not (R/r['path']).exists() for r in p['records']),gross_unique_allocated_bytes_retired=freed,compressed_allocated_bytes=objects['compressed_allocated_bytes'],filesystem_free_bytes_before=free_before,filesystem_free_bytes_after=shutil.disk_usage(R).free,same_user_processes_checked=processes,unique_data_discarded=False,all_archive_objects_verified=True)
 write(path.parent/('retirement-'+str(time.time_ns())+'.json'),summary);print(json.dumps(summary))

def restore(path,name):
 p=check_plan(path);rows=[r for r in p['records'] if r['path']==name];assert len(rows)==1;r=rows[0];file=R/name
 assert file.is_relative_to(R/'attestations') and not file.exists() and not file.is_symlink()
 objects=load(path.parent/'objects.json');blob=R/objects['objects'][r['sha256']]['path'];assert digest(blob)==objects['objects'][r['sha256']]['sha256']
 file.parent.mkdir(parents=True,exist_ok=True);assert file.parent.resolve()==file.parent
 temp=file.with_name('.sirius-restore-'+str(os.getpid()));assert not temp.exists()
 try:
  observed=hashlib.sha256();size=0
  with gzip.open(blob,'rb') as incoming,temp.open('xb') as stream:
   while data:=incoming.read(1048576):observed.update(data);size+=len(data);stream.write(data)
   stream.flush();os.fsync(stream.fileno())
  assert observed.hexdigest()==r['sha256'] and size==r['before']['bytes']
  os.chmod(temp,r['before']['mode']);os.utime(temp,ns=(r['before']['mtime_ns'],r['before']['mtime_ns']))
  assert meta(metadata(temp))==meta(r['before'])
  os.link(temp,file,follow_symlinks=False)
 finally:
  if temp.exists():temp.unlink()
 assert digest(file)==r['sha256'] and meta(metadata(file))==meta(r['before']);print(json.dumps(dict(restored=name,sha256=r['sha256'])))

parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['plan','pack','retire','restore']);parser.add_argument('plan',type=Path);parser.add_argument('--path')
a=parser.parse_args();path=a.plan.resolve();assert path.is_relative_to(R/'attestations/cleanup')
if a.mode=='plan':make_plan(path)
elif a.mode=='pack':pack(path)
elif a.mode=='retire':retire(path)
else:assert a.path;restore(path,a.path)
