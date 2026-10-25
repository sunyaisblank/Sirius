"""Preserve the sole original AN adoption-gate outcomes, without rerunning builds."""
from pathlib import Path
import hashlib,json,stat,subprocess,zipfile
R=Path.cwd();W=R/'out/native-cpu-local-identities';C=R/'out/native-cpu-local-identities-cohort';D=R/'out/native-cpu-local-identities-compatibility';J=lambda p:json.loads(p.read_text())
def identity(p):
 a=p.stat()
 with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 b=p.stat();assert (a.st_size,a.st_mode,a.st_mtime_ns)==(b.st_size,b.st_mode,b.st_mtime_ns)
 return dict(bytes=b.st_size,sha256=digest,mode=stat.S_IMODE(b.st_mode),mtime_ns=b.st_mtime_ns)
assert not subprocess.check_output(['git','status','--porcelain']).strip()
subprocess.run(['git','diff','--quiet','0109ab462edade524667d20303507aa09adcbaf6','HEAD','--','src','tests','cmake','CMakeLists.txt','CMakePresets.json'],check=True)
base=R/'attestations/qualification/86843614/native-cpu-local-identities-finite-compatibility/MANIFEST.json';prior=J(base)
assert prior['finite_component_passed'] and prior['original_C8_passed'] and prior['new_products']==9
recovery={};retired=J(base.parent/'RETIRED-DIAGNOSTIC-PRODUCTS.json');assert retired['status']=='complete';retired_paths={x['path'] for x in retired['members']}
for p in [base,*(base.parent/a['file'] for a in prior['archives']),base.parent/'RETIRED-DIAGNOSTIC-PRODUCTS.json',base.parent/'RETIRED-EXPANDED-LEDGERS.json']:
 assert subprocess.check_output(['git','show','@{upstream}:'+str(p.relative_to(R))])==p.read_bytes();recovery[str(p.relative_to(R))]=identity(p)
a,=[a for a in prior['archives'] if a['file']=='actual-products.zip'];current_products={}
with zipfile.ZipFile(base.parent/a['file']) as z:
 assert z.testzip() is None and set(z.namelist())==set(a['members'])
 for key,v in a['members'].items():
  data=z.read(key);assert len(data)==v['bytes'] and hashlib.sha256(data).hexdigest()==v['sha256'] and stat.S_IMODE(z.getinfo(key).external_attr>>16)==v['mode']
  p=Path(v['source']);relative=str(p.relative_to(R))
  if relative in retired_paths:assert not p.exists();continue
  assert identity(p)=={k:v[k] for k in ('bytes','sha256','mode','mtime_ns')} and p.read_bytes()==data;current_products[relative]=identity(p)
assert len(current_products)==6
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
  for row in o['observed_births']:
   if boot==o['boot_id']:
    try:t=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):t=None
    assert t!=str(row['saved_start_ticks'])
  owners[str(path.relative_to(R))]={k:o.get(k) for k in ('source','status','passed','child_exit','wall_seconds','peak_sampled_rss_kib','stop_reason')}
if all22 and all8 and frame_owner:assert len(owners)==20
source=J(C/'benefit/runtime/a1/owner.json')['source']['revision']
frame_passed=bool(frame_owner and frame_owner['passed'] and frame.get('accepted'))
scope='accepted-original-frame' if frame_passed else 'frame-timeout' if frame_owner and frame_owner['stop_reason']=='diagnostic_timeout' else 'adoption-gates-rejected' if not all22 else 'adoption-finite-outcome'
B=R/'attestations/qualification'/source[:8]/('native-cpu-local-identities-'+scope);assert not B.exists();B.mkdir(parents=True)
sources={str(p.relative_to(R)):p for p in C.rglob('*') if p.is_file()}
p=Path(__file__).resolve();sources[str(p.relative_to(R))]=p
path=B/'sole-adoption-cohort.zip';members={}
with zipfile.ZipFile(path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for key,p in sorted(sources.items()):
  before=identity(p);z.write(p,'files/'+key);assert identity(p)==before;members['files/'+key]=dict(source=str(p),**before)
with zipfile.ZipFile(path) as z:
 assert z.testzip() is None and set(z.namelist())==set(members)
 for key,v in members.items():
  data=z.read(key);assert len(data)==v['bytes'] and hashlib.sha256(data).hexdigest()==v['sha256'] and stat.S_IMODE(z.getinfo(key).external_attr>>16)==v['mode']
  assert Path(v['source']).read_bytes()==data and identity(Path(v['source']))=={k:v[k] for k in ('bytes','sha256','mode','mtime_ns')}
manifest=dict(source=source,accepted=frame_passed,original_frame_passed=frame_passed,all22_original_gates_pass=all22,all8_original_controls_pass=all8,original_frame_attempted=bool(frame_owner),full_scientific_and_release_qualification=False,owners=owners,new_compilations=0,new_links=0,new_products=0,reused_current_products=current_products,reused_pushed_recovery=recovery,archives=[dict(file=path.name,identity=identity(path),members=members)],
 scope='One first original adoption trial of AN after its previous finite component admission. Exactly original A1/B1/B2/A2 with three internal repeats, registered stock baseline and same H host, all22 thresholds and independent Decimal. Remaining7 only after all22PASS; original128x128/1SPP Kerr300s4GiB2GiB cold frame only after fresh all8/all22. No resampling, pooling/subtraction, warmup, changed counts/cache/governors/precision/settings/limits or observer in frame. Earlier AM failed frame is not rerun, and historical partial tile counts cannot establish causal speed/completion. AN is a distinct current implementation moving the exact complete-fact integer identity/refusal stores before the original C++ ABI fallback, with all writer/math/tuple/source/physical contracts unchanged. Original source/build/fixture/full38/C8 recovery and six current products verified against prior pushed archive; no compiler/link/workload is repeated here. Earlier AN component did not attempt these adoption gates; current actual outcomes alone are admitted above. Public source remains0109ab46. No full science/Mandatory/native-window/installed-release/exact-main readiness.')
if frame_owner:manifest['frame_terminal']={k:frame_owner.get(k) for k in ('source','passed','wall_seconds','peak_sampled_rss_kib','minimum_sampled_available_kib','stop_reason','child_exit')}
(B/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(destination=str(B.relative_to(R)),terminal_owners=len(owners),all22=all22,all8=all8,original_frame_passed=frame_passed,archive_bytes=manifest['archives'][0]['identity']['bytes'],archive_sha256=manifest['archives'][0]['identity']['sha256'])))
