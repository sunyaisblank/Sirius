"""Preserve the corrected-tag cohort and every original terminal outcome once."""
from pathlib import Path
import hashlib,json,stat,subprocess,zipfile
R=Path.cwd();W=R/'out/native-cpu-reachable-tag-writers';C=R/'out/native-cpu-reachable-tag-writers-cohort';D=R/'out/native-cpu-reachable-tag-writers-compatibility';J=lambda p:json.loads(p.read_text())
def identity(p):
 a=p.stat()
 with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 b=p.stat();assert (a.st_size,a.st_mode,a.st_mtime_ns)==(b.st_size,b.st_mode,b.st_mtime_ns)
 return dict(bytes=b.st_size,sha256=digest,mode=stat.S_IMODE(b.st_mode),mtime_ns=b.st_mtime_ns)
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert J(W/'fixture/result.json')['old_stale_tag_arithmetic_counterexamples']==12
assert J(W/'fixture/result.json')['fixed_writer_calls']==190512
assert J(W/'fixture/build-admission.json')['passed']
assert J(C/'control-8-acceptance.json')['passed'] and J(D/'result.json')['passed']
benefit=J(C/'benefit/readers/benefit.json') if (C/'benefit/readers/benefit.json').exists() else None
decimal=J(C/'benefit/readers/decimal-agreement.json') if (C/'benefit/readers/decimal-agreement.json').exists() else None
controls=J(C/'controls/readers/controls-1-8-acceptance.json') if (C/'controls/readers/controls-1-8-acceptance.json').exists() else None
frame=J(C/'cold-frame/result.json') if (C/'cold-frame/result.json').exists() else None
frame_owner=J(C/'cold-frame/frame/owner.json') if (C/'cold-frame/frame/owner.json').exists() else None
all22=bool(benefit and benefit['audit_passed'] and benefit['all22_numeric_pass'] and decimal and decimal['agrees_on_all22'] and decimal['all22_pass'])
all8=bool(controls and controls['all8_original_controls_pass'])
if frame_owner:assert frame and all22 and all8
owners={};boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
for root in (W,C,D):
 for path in sorted(root.rglob('owner.json')):
  o=J(path);assert o['status']=='terminal' and o['source_unchanged'] and o['whole_inputs_unchanged']
  assert not o['cleanup_errors'] and not o['remaining_owned_processes'] and all(x['absent'] for x in o['observed_births'])
  if path==W/'fixture/build/owner.json':
   a=J(W/'fixture/build-admission.json');assert not o['passed'] and a['original_owner_passed'] is False
   assert {k:identity(path)[k] for k in ('bytes','sha256')}==a['original_owner']
  for row in o['observed_births']:
   if boot==o['boot_id']:
    try:t=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):t=None
    assert t!=str(row['saved_start_ticks'])
  owners[str(path.relative_to(R))]={k:o.get(k) for k in ('source','status','passed','child_exit','wall_seconds','peak_sampled_rss_kib','stop_reason')}
source=J(C/'control-8/owner.json')['source']['revision']
frame_passed=bool(frame_owner and frame_owner['passed'] and frame.get('accepted'))
scope='accepted-original-frame' if frame_passed else 'frame-timeout' if frame_owner and frame_owner['stop_reason']=='diagnostic_timeout' else 'gated-finite-outcome'
B=R/'attestations/qualification'/source[:8]/('native-cpu-reachable-tag-writers-'+scope);assert not B.exists();B.mkdir(parents=True)
if all22 and all8 and frame_owner:assert len(owners)==20
compiled=J(W/'ispc-dense-compile-result.json');provider=J(W/'provider/provider-build-result.json');fixture=J(W/'fixture/build-result.json');diagnostic=J(D/'link-result.json')
products=[R/x['object']['path'] for x in compiled['objects']]+[R/x['product']['path'] for x in fixture['commands']]
products+=[R/provider[k] for k in ('library','ICD')];products+=[R/diagnostic['executable']['path']]
assert len(products)==len(set(products))==9
sources={str(p.relative_to(R)):p for root in (W,C,D) for p in root.rglob('*') if p.is_file()}
for row in J(W/'ispc-derivation.json')['stages'].values():
 p=R/row['derived_path'];sources[str(p.relative_to(R))]=p
p=Path(__file__).resolve();sources[str(p.relative_to(R))]=p
def archive(name,files):
 path=B/name;members={}
 with zipfile.ZipFile(path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for key,p in sorted(files.items()):
   before=identity(p);z.write(p,'files/'+key);assert identity(p)==before;members['files/'+key]=dict(source=str(p),**before)
 with zipfile.ZipFile(path) as z:
  assert z.testzip() is None and set(z.namelist())==set(members)
  for key,v in members.items():
   data=z.read(key);assert len(data)==v['bytes'] and hashlib.sha256(data).hexdigest()==v['sha256']
   assert stat.S_IMODE(z.getinfo(key).external_attr>>16)==v['mode']
   assert Path(v['source']).read_bytes()==data and identity(Path(v['source']))=={k:v[k] for k in ('bytes','sha256','mode','mtime_ns')}
 return dict(file=path.name,identity=identity(path),members=members)
archives=[archive('source-build-and-trial.zip',sources),archive('actual-products.zip',{str(p.relative_to(R)):p for p in products})]
references={}
for name in ['attestations/qualification/0ce77927/native-cpu-normalization-tags-frame-timeout/MANIFEST.json','attestations/qualification/5bc67064/native-cpu-predicated-refinement-finite-compatibility/MANIFEST.json','attestations/qualification/f878a18f/native-cpu-normalized-arithmetic-frame-timeout/MANIFEST.json']:
 p=R/name;assert subprocess.check_output(['git','show','@{upstream}:'+name])==p.read_bytes();references[name]=identity(p)
manifest=dict(source=source,accepted=frame_passed,original_frame_passed=frame_passed,all22_original_gates_pass=all22,all8_original_controls_pass=all8,
 original_frame_attempted=bool(frame_owner),full_scientific_and_release_qualification=False,owners=owners,archives=archives,new_products=9,
 production_raw_compilations=3,production_CPP_compilations=0,diagnostic_ISPC_compilations=1,diagnostic_CPP_compilations=1,reused_provenance=references,
 superseded_claim='AJ and AK raw tag invalidation/refresh followed unconditional returns. Their coherent-tag claims are disproved by twelve actual-body stale-tag arithmetic counterexamples. Prior finite numerical observations remain historical; no coherent-tag admission is inherited from them. Only exact AK arithmetic/adapter objects and separately pinned unaffected components are reused.',
 inherited_template_labels='Some original remaining-control scope strings retain P host/Zprovider wording. Actual sealed action executable/ICD/provider identities, observed mappings and original readers bind the reused H host and corrected AL candidate; A1/A2 separately bind the registered stock baseline. The inherited prose is not identity authority. No original record is rewritten.',
 scope='AL moves twelve metadata updates before their original returns; exact inverse restores AJ and canonical9380tuples/public payload/math/guards/masks/order remain. Actual writer regression190512calls/22861440bank/4572288tag words plus12 old arithmetic counterexamples passes. Original fixture build-reader failure remains failed with explicit unchanged-link/GOT/relaxation-target admission. Current provider894function/14target joins, full38/187960words and original C8/all32 pass. Any original ABBA22/Decimal/all8/frame results are recorded above only if actually completed. No full science/Mandatory/native-window/installed-release/public-main readiness claim.')
if frame_owner:manifest['frame_terminal']={k:frame_owner.get(k) for k in ('source','passed','wall_seconds','peak_sampled_rss_kib','minimum_sampled_available_kib','stop_reason','child_exit')}
(B/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(destination=str(B.relative_to(R)),terminal_owners=len(owners),new_products=9,all22=all22,all8=all8,original_frame_passed=frame_passed,archives=[dict(file=a['file'],bytes=a['identity']['bytes'],sha256=a['identity']['sha256']) for a in archives])))
