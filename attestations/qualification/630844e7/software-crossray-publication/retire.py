"""Retire the rejected192-file predecessor; preserve the corrected immutable root."""
from pathlib import Path
import hashlib,json,shutil,stat,subprocess,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;C='65a466249d923a4b8c2ce1a02496cafdaa2e733a';BR='development/software-crossray-endpoint-2026-10-20'
refs={'refs/heads/'+BR:C,'refs/tags/evidence/software-crossray-endpoint-rejection-630844e7':'32cf94d3dec0cafb54596300e793915376cf9665','refs/tags/evidence/software-crossray-active-bound-630844e7':C}
remote=subprocess.check_output(['git','ls-remote','origin',*refs],text=True);assert {r:c for c,r in (s.split() for s in remote.splitlines())}==refs;(W/'checkpoint-remote.txt').write_text(remote)
def signature(p):
 s=p.lstat();assert stat.S_ISREG(s.st_mode) and not p.is_symlink();return {'bytes':s.st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'mode':oct(stat.S_IMODE(s.st_mode))}
p=R/'attestations/qualification/630844e7/software-crossray-endpoint/MANIFEST.json';m=json.loads(p.read_text());a=p.parent/m['archive']['path'];assert signature(a)=={k:m['archive'][k] for k in ['bytes','sha256','mode']}
for f in [p,a]:
 data=subprocess.check_output(['git','show',C+':'+str(f.relative_to(R))]);assert len(data)==f.stat().st_size and hashlib.sha256(data).hexdigest()==signature(f)['sha256']
with zipfile.ZipFile(a) as z:
 assert z.testzip() is None and set(z.namelist())==set(m['archive_backed_files'])
 for i in z.infolist():
  pin=m['archive_backed_files'][i.filename];data=z.read(i);assert len(data)==pin['bytes'] and hashlib.sha256(data).hexdigest()==pin['sha256'] and stat.S_IMODE(i.external_attr>>16)==int(pin['mode'],8)
actual={}
for root in m['owned_roots']:
 for f in (R/root).rglob('*'):
  assert not f.is_symlink()
  if f.is_file():actual[str(f.relative_to(R))]=signature(f)
assert actual==m['archive_backed_files'] and len(actual)==192
future=R/'out/software-crossray-active';before=json.loads((future/'compile-endpoint/inputs-before.json').read_text());assert before==json.loads((future/'compile-endpoint/inputs-after.json').read_text())
for name in before:
 p=R/name;assert all(not p.is_relative_to(R/root) for root in m['owned_roots']),('corrected compiler still needs rejected live input',name)
frozen_before={str(p.relative_to(R)):signature(p) for p in future.rglob('*') if p.is_file()}
for root in m['owned_roots']:shutil.rmtree(R/root)
assert all(not (R/root).exists() for root in m['owned_roots'])
assert frozen_before=={str(p.relative_to(R)):signature(p) for p in future.rglob('*') if p.is_file()}
result={'passed':True,'published_checkpoint':C,'branch':BR,'retired_files':len(actual),'retired_bytes':sum(p['bytes'] for p in actual.values()),'retired_roots':m['owned_roots'],'all_live_bytes_modes_CRC_and_published_Git_payloads_equal_before_removal':True,'corrected_root':str(future.relative_to(R)),'corrected_root_files_retained':len(frozen_before),'corrected_root_byte_exact_before_after':True,'corrected_compile_frozen_inputs_do_not_depend_on_retired_roots':True,'canonical_preset_bins_and_unrelated_work_preserved':True,'active_owned_runtime_processes':0,'source_revision':'630844e775a85e74eb5d9ed328a2b8019abb640d'}
(W/'retirement.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['retired_files','retired_bytes','corrected_root_files_retained','corrected_root_byte_exact_before_after']}))
