"""Recover the finite AK source/products without claiming benefit or a frame."""
from pathlib import Path
import hashlib,json,os,stat,subprocess,zipfile
R=Path.cwd();W=R/'out/native-cpu-predicated-refinement';C=R/'out/native-cpu-predicated-refinement-cohort';D=R/'out/native-cpu-predicated-refinement-compatibility'
J=lambda p:json.loads(p.read_text())
def identity(p):
 a=p.stat()
 with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 b=p.stat();assert (a.st_size,a.st_mtime_ns,a.st_mode)==(b.st_size,b.st_mtime_ns,b.st_mode)
 return dict(bytes=b.st_size,sha256=digest,mode=stat.S_IMODE(b.st_mode),mtime_ns=b.st_mtime_ns)
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert J(C/'control-8-acceptance.json')['passed'] and J(C/'affected-control/result.json')['passed'] and J(D/'result.json')['passed']
assert J(W/'cpp/run/stdout.log')['compared_words']==44452800 and J(W/'cpp/run/stdout.log')['tag_words']==8890560
assert not (C/'cold-frame').exists() and not (C/'benefit').exists()
owners={};boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
for root in (W,C,D):
 for p in sorted(root.rglob('owner.json')):
  o=J(p);assert o['status']=='terminal' and o['passed'] and o['source_unchanged'] and o['whole_inputs_unchanged']
  assert not o['cleanup_errors'] and not o['remaining_owned_processes']
  assert o['observed_births'] and all(x['absent'] for x in o['observed_births'])
  for row in o['observed_births']:
   if boot==o['boot_id']:
    try:current=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):current=None
    assert current!=str(row['saved_start_ticks'])
  owners[str(p.relative_to(R))]={k:o.get(k) for k in ('source','status','passed','child_exit','wall_seconds','peak_sampled_rss_kib','stop_reason')}
assert len(owners)==8
source=J(C/'control-8/owner.json')['source']['revision']
B=R/'attestations/qualification'/source[:8]/'native-cpu-predicated-refinement-finite-compatibility';assert not B.exists();B.mkdir(parents=True)
cpp=J(W/'cpp-build-result.json');provider=J(W/'provider/provider-build-result.json');linked=J(D/'link-result.json')
assert cpp['new_CPP_compilations']==1 and cpp['new_relocatable_links']==1
products=[R/x['object']['path'] for x in cpp['compilations']]
products+=[R/cpp['merged_adapter']['path'],R/cpp['executable']['path'],*(R/provider[k] for k in ('library','ICD')),R/linked['executable']['path']]
assert len(products)==len(set(products))==6
sources={str(p.relative_to(R)):p for root in (W,C,D) for p in root.rglob('*') if p.is_file()}
sources[str(Path(__file__).resolve().relative_to(R))]=Path(__file__).resolve()
def archive(name,files):
 path=B/name;members={}
 with zipfile.ZipFile(path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for name,p in sorted(files.items()):
   value=identity(p);z.write(p,'files/'+name);assert identity(p)==value
   members['files/'+name]=dict(source=str(p),**value)
 with zipfile.ZipFile(path) as z:
  assert z.testzip() is None and set(z.namelist())==set(members)
  for name,value in members.items():
   data=z.read(name);assert len(data)==value['bytes'] and hashlib.sha256(data).hexdigest()==value['sha256']
   assert stat.S_IMODE(z.getinfo(name).external_attr>>16)==value['mode']
   assert (R/Path(value['source'])).read_bytes()==data
   assert identity(Path(value['source']))=={k:value[k] for k in ('bytes','sha256','mode','mtime_ns')}
 return dict(file=path.name,identity=identity(path),members=members)
archives=[archive('source-build-and-controls.zip',sources),archive('actual-products.zip',{str(p.relative_to(R)):p for p in products})]
references={}
for name in ['attestations/qualification/0ce77927/native-cpu-normalization-tags-frame-timeout/MANIFEST.json','attestations/qualification/f878a18f/native-cpu-normalized-arithmetic-frame-timeout/MANIFEST.json','attestations/qualification/db09db3f/native-cpu-steady-profile/MANIFEST.json']:
 p=R/name;assert subprocess.check_output(['git','show','@{upstream}:'+name])==p.read_bytes();references[name]=identity(p)
manifest=dict(source=source,accepted=False,finite_compatibility_passed=True,frame_executed=False,benefit_ABBA_executed=False,
 full_scientific_and_release_qualification=False,new_CPP_compilations=1,reused_raw_objects=3,new_products=6,
 owners=owners,archives=archives,reused_provenance=references,
 scope='AK six empty-work terminal returns in C++ Divide/Sqrt; source inverse restores AJ exactly. All nonempty refinements retain original formulas/order/three correction steps. One helper compilation, original driver/raw3/scalar3/pure admission/adapter/H host reused. Original370440call/44452800payload-word/8890560tag-word fixture, actual892function/14target provider joins, fresh full38/187960word comparison, original C8/all32 and affected original connected-interval C5 pass. No current ABBA/all8/frame performance acceptance, public adoption, science/Mandatory/native-window/installed-release claim. Historical Z all8/all22 and earlier rejected frames retain their original scopes. Further300s frame not justified by this small component change alone.')
(B/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(destination=str(B.relative_to(R)),terminal_owners=len(owners),products=len(products),archives=[dict(file=x['file'],bytes=x['identity']['bytes'],sha256=x['identity']['sha256']) for x in archives])))
