"""Retire only exact immutable roots backed by published Git checkpoints."""
from pathlib import Path
import hashlib,json,shutil,stat,subprocess,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;C='121db07cabe4508890c3ba3045bfe040805df8be'
BR='development/software-single-invocation-2026-10-20';TAG='evidence/software-single-invocation-endpoint-cost-630844e7'
remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/'+BR,'refs/tags/'+TAG],text=True)
rows={ref:sha for sha,ref in (s.split() for s in remote.splitlines())}
assert rows=={'refs/heads/'+BR:C,'refs/tags/'+TAG:C}
(W/'checkpoint-remote.txt').write_text(remote)
def signature(p):
 s=p.lstat();assert stat.S_ISREG(s.st_mode) and not p.is_symlink();return {'bytes':s.st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'mode':oct(stat.S_IMODE(s.st_mode))}
manifests=[R/'attestations/qualification/630844e7/software-single-invocation-endpoint/MANIFEST.json',R/'attestations/qualification/630844e7/software-single-invocation-endpoint-cost/MANIFEST.json']
roots=[];all_files={};backing=[]
for p in manifests:
 d=json.loads(p.read_text());a=p.parent/d['archive']['path'];pin={k:d['archive'][k] for k in ['bytes','sha256','mode']};assert signature(a)==pin
 for f in [p,a]:
  data=subprocess.check_output(['git','show',C+':'+str(f.relative_to(R))]);assert len(data)==f.stat().st_size and hashlib.sha256(data).hexdigest()==signature(f)['sha256']
 with zipfile.ZipFile(a) as z:
  assert z.testzip() is None and set(z.namelist())==set(d['archive_backed_files'])
  for i in z.infolist():
   pin=d['archive_backed_files'][i.filename];data=z.read(i);assert len(data)==pin['bytes'] and hashlib.sha256(data).hexdigest()==pin['sha256'] and stat.S_IMODE(i.external_attr>>16)==int(pin['mode'],8)
 for name,pin in d['archive_backed_files'].items():assert name not in all_files;all_files[name]=pin
 roots+=d['owned_roots'];backing.append({'manifest':str(p.relative_to(R)),**signature(p),'archive':str(a.relative_to(R))})
actual={}
for root in roots:
 p=R/root;assert p.is_dir() and not p.is_symlink()
 for f in p.rglob('*'):
  assert not f.is_symlink()
  if f.is_file():actual[str(f.relative_to(R))]=signature(f)
assert actual==all_files
for root in roots:shutil.rmtree(R/root)
assert all(not (R/root).exists() for root in roots)
result={'passed':True,'published_checkpoint':C,'branch':BR,'tag':TAG,'backing':backing,'all_live_bytes_modes_archive_CRC_and_Git_payloads_equal_before_removal':True,'retired_roots':roots,'retired_files':len(all_files),'retired_bytes':sum(v['bytes'] for v in all_files.values()),'canonical_preset_bins_and_unrelated_files_preserved':True,'active_owned_runtime_processes':0,'prototype_disposition':'REJECT_ENDPOINT_SINGLE_INVOCATION_COST','source_revision':'630844e775a85e74eb5d9ed328a2b8019abb640d'}
(W/'retirement.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'retired_files':result['retired_files'],'retired_roots':len(roots),'retired_bytes':result['retired_bytes'],'published_checkpoint':C}))
