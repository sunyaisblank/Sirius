from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
scratch=ROOT/'out/native-homogeneous-opcode-lowering';archive=ROOT/'attestations/native-vulkan/f02ac83/homogeneous-opcode-lowering'
receipt=ROOT/'attestations/cleanup/2026-10-10/homogeneous-opcode-native-query-finished.json'
assert not receipt.exists() and scratch.is_dir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def document(p):return json.loads(p.read_text(encoding='utf-8-sig'))
manifest=archive/'manifest.json';assert sha(manifest)=='b60c39cb82ebf212c84ecb5c66c39eb81e5360fdf05fc60b0c5033933b87b977'
m=document(manifest);records=m['files'];assert len(records)==33 and sum(q['bytes'] for q in records)==11118945
for q in records:
 for root in (archive,scratch):
  p=root/q['path'];assert p.is_file() and not p.is_symlink() and p.stat().st_size==q['bytes'] and sha(p)==q['sha256'],str(p)
expected={q['path'] for q in records};extra='cleanup-static-query.json'
actual={str(p.relative_to(scratch)) for p in scratch.rglob('*') if p.is_file()}
assert actual==expected|{extra},actual-expected
assert not any(p.is_symlink() for p in scratch.rglob('*'))
head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();assert head=='bd82e523227a1b7223f43984af10291711925e44'
assert not subprocess.check_output(['git','status','--porcelain'])
assert subprocess.check_output(['git','rev-list','--left-right','--count','HEAD...@{upstream}'],text=True).strip()=='0\t0'
assert subprocess.check_output(['git','show',head+':src/sirius/kernels/retained_transport.slang'])==(archive/'candidate-source.slang').read_bytes()
assert subprocess.check_output(['git','ls-remote','origin','refs/heads/development/retained-dense-step-control-2026-10-06'],text=True).split()[0]==head
publication=document(WORK/'static-query-publication.json');assert publication['readback']['body']==(WORK/'static-query-publication.txt').read_text()
windows=document(scratch/extra);bridge=document(WORK/'cleanup-static-query-bridge.json')
assert bridge['accepted'] and bridge['returncode']==0 and not bridge['outer_stop'] and not bridge['errors'] and bridge['linux_child_birth_absent']
assert windows['owned_births_absent'] and windows['all_two_module_seals_match'] and not windows['consumer_scan']['matches'] and not windows['consumer_scan']['accessible_compilers']
def stat(pid):return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()
# Exclude only this read-only inspector and its exact currently observed ancestors.
ancestors=[];current=os.getpid()
while current>0:
 try:f=stat(current)
 except (FileNotFoundError,ProcessLookupError):break
 ancestors.append({'pid':current,'start_ticks':f[19]});current=int(f[1])
excluded={(q['pid'],q['start_ticks']) for q in ancestors}
refs=[str(scratch),'out/native-homogeneous-opcode-lowering'];matches=[];unreadable={};processes=0
for folder in Path('/proc').iterdir():
 if not folder.name.isdecimal():continue
 pid=int(folder.name)
 try:f=stat(pid)
 except (FileNotFoundError,ProcessLookupError):continue
 except OSError:unreadable['stat']=unreadable.get('stat',0)+1;continue
 if (pid,f[19]) in excluded:continue
 processes+=1
 fields=[]
 for field in ('cmdline','exe','cwd'):
  try:
   value=(folder/field).read_bytes().replace(b'\0',b' ').decode(errors='replace') if field=='cmdline' else os.readlink(folder/field)
   if any(word in value for word in refs):fields.append(field)
  except (FileNotFoundError,ProcessLookupError):pass
  except OSError:unreadable[field]=unreadable.get(field,0)+1
 try:
  for fd in (folder/'fd').iterdir():
   try:
    value=os.readlink(fd)
    if any(word in value for word in refs):fields.append('fd:'+fd.name)
   except (FileNotFoundError,ProcessLookupError):pass
   except OSError:unreadable['fd_target']=unreadable.get('fd_target',0)+1
 except (FileNotFoundError,ProcessLookupError):pass
 except OSError:unreadable['fd_directory']=unreadable.get('fd_directory',0)+1
 if fields:matches.append({'pid':pid,'start_ticks':f[19],'fields':fields})
assert not matches,matches
old=document(scratch/'native-bridge.json');pid=old['linux_owner_pid']
try:birth=stat(pid)[19]
except (FileNotFoundError,ProcessLookupError):birth=None
assert birth!=old['linux_owner_start_ticks']
linux={'processes':processes,'unreadable_field_counts':unreadable,'matches':matches,'excluded_exact_inspector_ancestors':ancestors,'saved_query_bridge':{'pid':pid,'saved_start_ticks':old['linux_owner_start_ticks'],'current_start_ticks':birth,'owned_birth_absent':True},'scope':'Accessible Linux cmdline/executable/cwd/fd references and exact saved bridge birth; no hidden/unreadable/global consumer absence guarantee.'}
# Reverify the frozen payloads immediately before root-only removal.
for q in records:
 for root in (archive,scratch):
  p=root/q['path'];assert p.stat().st_size==q['bytes'] and sha(p)==q['sha256']
proof=receipt.with_suffix('').with_name('homogeneous-opcode-native-query-final-proofs');assert not proof.exists();proof.mkdir()
for name,p in [('windows-terminal.json',scratch/extra),('windows-terminal-bridge.json',WORK/'cleanup-static-query-bridge.json'),('query-final-observer.py',WORK/'run_query_cleanup_audit.py'),('cleanup.py',Path(__file__))]:
 shutil.copyfile(p,proof/name);assert sha(p)==sha(proof/name)
(proof/'linux-consumer-audit.json').write_text(json.dumps(linux,indent=2)+'\n')
result={'completed_utc':datetime.now(timezone.utc).isoformat(),'source_revision':head,'committed_and_pushed':True,'branch':'development/retained-dense-step-control-2026-10-06','tracking':'0/0','archive':{'path':str(archive.relative_to(ROOT)),'files':33,'bytes':11118945,'manifest_sha256':sha(manifest),'whole_payloads_reverified':True},'independent_archive_review':'/root/cleanup_review rehashed all33 frozen files and confirmed archived source exactly in pushedbd82e523; no preservation blocker','static_disposition_comment':publication['readback']['html_url'],'linux':linux,'windows':windows,'final_proofs':str(proof.relative_to(ROOT)),'removal_target':str(scratch.relative_to(ROOT)),'finished_scratch_removed':False,'active_trial_preserved':'out/transport-homogeneous-opcode-trial','scope':'Completed zero-dispatch native compiler query only. Useful source and concise findings pushed/published; raw evidence local/ignored without remote-backup guarantee. Candidate controls/retention/adoption and original full frames/full qualification remain pending.'}
receipt.write_text(json.dumps(result,indent=2)+'\n')
shutil.rmtree(scratch);assert not scratch.exists();result['finished_scratch_removed']=True;receipt.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'removed':str(scratch.relative_to(ROOT)),'archive_files':33,'archive_bytes':11118945,'linux_processes':processes,'linux_unreadable':unreadable,'windows_processes':windows['consumer_scan']['processes'],'windows_hidden':windows['consumer_scan']['hidden_field_records'],'active_trial_preserved':True}))
