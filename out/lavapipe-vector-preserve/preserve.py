"""Preserve the rejected private compiler candidate and closed original trial."""
from pathlib import Path
import hashlib,importlib.util,json,os,shutil,stat,subprocess,tarfile,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;D=R/'attestations/qualification/8c2be67d/lavapipe-vector-workgroup-rejected'
assert not D.exists();D.mkdir(parents=True)
OUT=['out/lavapipe-vector-'+n for n in ['workgroup','producer','candidate','inspect-source','link','inspect','capture','control-source','validation-build','runtime','readers','preserve']]
BIN=['bin/linux-gcc/lavapipe-vector-'+n for n in ['staged','inspect','validation']]
BUILD=R/'bin/linux-gcc/lavapipe-vector-workgroup'
def digest(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def identity(p):return {'bytes':p.stat().st_size,'sha256':digest(p),'mode':stat.S_IMODE(p.stat().st_mode)}
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
v=json.loads((R/'out/lavapipe-vector-readers/benefit.json').read_text());assert v['audit_passed'] and not v['benefit_passed'] and (v['checks_passed'],v['checks_total'])==(16,22)
assert json.loads((R/'out/lavapipe-vector-readers/decimal-agreement.json').read_text())['agrees_on_all22']
subprocess.run(['git','diff','--quiet','6d2af336b12dc430792bd5a83c04fb5632868c04','--','src','scripts','tests','docs','CMakeLists.txt','CMakePresets.json'],check=True)
spec=importlib.util.spec_from_file_location('closed_owner',R/'out/lavapipe-vector-producer/owned_guard.py');guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
owners={};births={}
for root in OUT:
 for p in (R/root).rglob('owner.json'):
  if '/upstream/' in str(p):continue
  o=json.loads(p.read_text());assert o['passed'] and o['status']=='terminal' and o['remaining_owned_processes']==[] and not o['cleanup_errors']
  assert o['stop_reason'] is None and o['source_unchanged'] and o['whole_inputs_unchanged']
  owners[str(p.relative_to(R))]={'identity':identity(p),'wall_seconds':o['wall_seconds'],'peak_sampled_rss_kib':o['peak_sampled_rss_kib'],'source_revision':o['source']['revision']}
  for item in o['observed_births']:
   assert item['absent'];pid,ticks=item['pid'],item['saved_start_ticks']
   try:now=guard.proc(pid)['start_ticks']
   except (FileNotFoundError,ProcessLookupError):now=None
   assert now!=ticks;(births.setdefault(str(pid)+':'+str(ticks),{'pid':pid,'saved_start_ticks':ticks,'current_start_ticks':now,'absent':True}))
assert len(owners)==17,len(owners)
upstream=R/'out/lavapipe-vector-workgroup/upstream/mesa-25.2.8';archive=R/'out/lavapipe-vector-workgroup/mesa-25.2.8.tar.xz'
admission=json.loads((R/'out/lavapipe-vector-workgroup/upstream-admission.json').read_text());assert digest(archive)==admission['archive_sha256'] and archive.stat().st_size==admission['archive_bytes']
derived={};changed=set(admission['changed_files'])
with tarfile.open(archive,'r:xz') as source:
 for member in source:
  if not member.isfile():continue
  rel=Path(member.name).relative_to('mesa-25.2.8');p=upstream/rel
  raw=(R/'out/lavapipe-vector-workgroup/candidate'/rel).read_bytes() if str(rel) in changed else source.extractfile(member).read()
  assert p.read_bytes()==raw,str(rel);derived[str(p.relative_to(R))]=identity(p)
assert len(derived)==11100 and set(derived)=={str(p.relative_to(R)) for p in upstream.rglob('*') if p.is_file()}
shutil.copy2(archive,D/archive.name)
save(D/'upstream-reconstruction.json',{'archive':identity(D/archive.name),'official_URL':admission['url'],'official_checksum_source':admission['official_checksum_source'],'patched_files':sorted(changed),'all11100_actual_inputs_reconstruct_exactly':True,'files':derived})
files={};archives={}
def pack(name,paths):
 target=D/name
 with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
  for p in sorted(set(paths)):
   rel=str(p.relative_to(R));assert rel not in files
   assert p.is_file() and not p.is_symlink(),rel
   files[rel]={**identity(p),'archive':name};z.write(p,rel)
 with zipfile.ZipFile(target) as z:
  assert z.testzip() is None
  for info in z.infolist():
   old=files[info.filename];raw=z.read(info)
   assert len(raw)==old['bytes'] and hashlib.sha256(raw).hexdigest()==old['sha256']
   assert stat.S_IMODE(info.external_attr>>16)==old['mode']
 archives[name]={**identity(target),'members':len(z.infolist())}
diagnostics=[]
for root in OUT:
 for p in (R/root).rglob('*'):
  if p.is_file() and not p.is_relative_to(upstream) and p!=archive:diagnostics.append(p)
pack('diagnostics-and-trial.zip',diagnostics)
products=[p for root in BIN for p in (R/root).rglob('*') if p.is_file()]
products += [R/'bin/linux-gcc/tests/backend/portable_binary32_reference.bin']
for rel in ['build.ninja','compile_commands.json','meson-private/coredata.dat','meson-private/cmd_line.txt','meson-logs/meson-log.txt','src/gallium/drivers/llvmpipe/libllvmpipe.a.p/lp_state_cs.c.o']:
 p=BUILD/rel;assert p.is_file();products.append(p)
pack('providers-and-probes.zip',products)
save(D/'MANIFEST.json',{'schema':'sirius-private-lavapipe-vector-rejection-v1','canonical_Sirius_source_revision':'6d2af336b12dc430792bd5a83c04fb5632868c04',
 'private_provider_source_revision':'5df02b1175805c18003f09959908ed7c14ba5837','runtime_revision':'8c2be67d696269679821cea95595a54da4b9dafc',
 'source_and_helper_preservation_revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 'roots':OUT+BIN+[str(BUILD.relative_to(R))],'archives':archives,'official_source_tar':identity(D/archive.name),'files':files,
 'terminal_owners':owners,'current_births':births,'reconstructable_upstream_files':11100,
 'build_cache_scope':'Remaining private Meson object/dependency caches are disposable reproducible outputs; both actual providers, matched complete llvmpipe archives, selected baseline/candidate objects, actual observer/probe ELFs, generated SPVs, configure/link commands and closed producers are backed. System tools/libraries are frozen identities, not backups.',
 'decision':'REJECT16_OF22','whole_passed':3,'whole_total':4,'stage_passed':13,'stage_total':16,'API_passed':0,'API_total':2,
 'no_remaining7_controls_frame_retry_or_adoption':True,'main_and_PR62_source_unchanged':True,
 'scope':'Closed private compiler experiment only. Original arithmetic and C8/sole ABBA correctness passed; fixed original22 benefit failed. No full Mandatory/scientific/native/provider/frame/installed-release qualification.'})
print(json.dumps({'preserved':True,'owners':len(owners),'absent_births':len(births),'files':len(files),'upstream_files_reconstructable':len(derived),'archives':archives,'official_tar':identity(D/archive.name)},indent=2))
