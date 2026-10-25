"""Preserve complete value facts, identity stores and each original terminal outcome once."""
from pathlib import Path
import hashlib,json,stat,subprocess,zipfile
R=Path.cwd();W=R/'out/native-cpu-value-fact-stores';C=R/'out/native-cpu-value-fact-stores-cohort';D=R/'out/native-cpu-value-fact-stores-compatibility';S=R/'out/native-cpu-value-facts';J=lambda p:json.loads(p.read_text())
def identity(p):
 a=p.stat()
 with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 b=p.stat();assert (a.st_size,a.st_mode,a.st_mtime_ns)==(b.st_size,b.st_mode,b.st_mtime_ns)
 return dict(bytes=b.st_size,sha256=digest,mode=stat.S_IMODE(b.st_mode),mtime_ns=b.st_mtime_ns)
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert J(W/'fixture/result.json')['old_stale_tag_arithmetic_counterexamples']==24
assert J(W/'fixture/result.json')['fixed_writer_calls']==254016
assert J(W/'fixture/build/owner.json')['passed'] and J(W/'fixture/linked-code.json')['passed']
assert J(C/'control-8-acceptance.json')['passed'] and J(D/'result.json')['passed']
for root in (S,W):
 result=J(root/'cpp/run/stdout.log');assert result['passed'] and result['calls']==555660 and result['compared_words']==66679200 and result['tag_words']==13335840 and result['active_lane_results']==1666980 and result['case_tuples']==49
 assert result['input_alias_and_inactive_words_checked'] and result['FP_controls_unchanged'] and J(root/'cpp-code-inspection.json')['passed']
benefit=J(C/'benefit/readers/benefit.json') if (C/'benefit/readers/benefit.json').exists() else None
decimal=J(C/'benefit/readers/decimal-agreement.json') if (C/'benefit/readers/decimal-agreement.json').exists() else None
controls=J(C/'controls/readers/controls-1-8-acceptance.json') if (C/'controls/readers/controls-1-8-acceptance.json').exists() else None
frame=J(C/'cold-frame/result.json') if (C/'cold-frame/result.json').exists() else None
frame_owner=J(C/'cold-frame/frame/owner.json') if (C/'cold-frame/frame/owner.json').exists() else None
all22=bool(benefit and benefit['audit_passed'] and benefit['all22_numeric_pass'] and decimal and decimal['agrees_on_all22'] and decimal['all22_pass'])
all8=bool(controls and controls['all8_original_controls_pass'])
if frame_owner:assert frame and all22 and all8
owners={};boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
for root in (S,W,C,D):
 for path in sorted(root.rglob('owner.json')):
  o=J(path);assert o['status']=='terminal' and o['source_unchanged'] and o['whole_inputs_unchanged']
  assert not o['cleanup_errors'] and not o['remaining_owned_processes'] and all(x['absent'] for x in o['observed_births'])
  assert o['passed'] or path==C/'cold-frame/frame/owner.json'
  for row in o['observed_births']:
   if boot==o['boot_id']:
    try:t=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):t=None
    assert t!=str(row['saved_start_ticks'])
  owners[str(path.relative_to(R))]={k:o.get(k) for k in ('source','status','passed','child_exit','wall_seconds','peak_sampled_rss_kib','stop_reason')}
source=J(C/'control-8/owner.json')['source']['revision']
frame_passed=bool(frame_owner and frame_owner['passed'] and frame.get('accepted'))
scope='accepted-original-frame' if frame_passed else 'frame-timeout' if frame_owner and frame_owner['stop_reason']=='diagnostic_timeout' else 'gated-finite-outcome'
B=R/'attestations/qualification'/source[:8]/('native-cpu-value-fact-stores-'+scope);assert not B.exists();B.mkdir(parents=True)
if all22 and all8 and frame_owner:assert len(owners)==24
compiled=J(W/'ispc-dense-compile-result.json');provider=J(W/'provider/provider-build-result.json');fixture=J(W/'fixture/build-result.json');diagnostic=J(D/'link-result.json')
products=[R/x['object']['path'] for x in compiled['objects']]+[R/x['product']['path'] for x in fixture['commands']]
products+=[R/provider[k] for k in ('library','ICD')];products+=[R/diagnostic['executable']['path']]
for root in (S,W):
 cpp=J(root/'cpp-build-result.json');assert cpp['passed']
 products += [R/x['object']['path'] for x in cpp['compilations']] + [R/cpp[k]['path'] for k in ('executable','merged_adapter')]
assert len(products)==len(set(products))==17
sources={str(p.relative_to(R)):p for root in (S,W,C,D) for p in root.rglob('*') if p.is_file()}
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
for name in ['attestations/qualification/893125b9/native-cpu-reachable-tag-writers-frame-timeout/MANIFEST.json','attestations/qualification/f878a18f/native-cpu-normalized-arithmetic-frame-timeout/MANIFEST.json']:
 p=R/name;assert subprocess.check_output(['git','show','@{upstream}:'+name])==p.read_bytes();references[name]=identity(p)
manifest=dict(source=source,accepted=frame_passed,original_frame_passed=frame_passed,all22_original_gates_pass=all22,all8_original_controls_pass=all8,
 original_frame_attempted=bool(frame_owner),full_scientific_and_release_qualification=False,owners=owners,archives=archives,new_products=17,
 production_raw_compilations=3,production_CPP_compilations=1,diagnostic_ISPC_compilations=1,diagnostic_CPP_compilations=2,
 supporting_component_CPP_compilations=2,reused_provenance=references,
 legacy_defect='AL corrected twelve unreachable AJ/AK writer metadata updates. AM extends that tested reachable protocol with complete exact current value facts. The24 actual old-writer witnesses comprise12 historical legacy normalization-bit witnesses and12 rich-schema/current-math extension witnesses; they are not24 historical AJ-binary errors.',
 inherited_template_labels='The writer fixture scope retains an AK math label, but its sealed linked object is actual AM probe/retained8.o SHA007dcf466adf1e513c6c60f7a67ee3692dfe7f2d6d13fd5602609f4edb69ff92. The provider source contract explicitly records this. Actual sealed products/code joins are authority; no original record is rewritten.',
 scope='AM caches exact stored zero/plus-one/minus-one/valid facts in the existing metadata word and uses original Add/Sub/Multiply identity/refusal precedence before unnecessary payload loads and normalization checks. Complete stores refresh; partial stores clear; unknown/legacy metadata conservatively falls back. Raw complete metadata guards identity stores; temporary local operand facts never authorize raw bank stores. Numeric signed zero and noncanonical nonzero valid values preserve original payload behavior. General arithmetic ordering/rounding/guards, Divide/Sqrt eligibility and three refinements, canonical9380tuples/public five-word payload, masks, alias snapshots, strict FP/ISA and scheduler stay unchanged. No intermediate FP exception-status equivalence is claimed. Supporting basic-facts and final identity-store pieces each pass555660calls/66679200bank/13335840metadata words. Actual writer fixture254016calls/30481920bank/6096384metadata words plus24 defined stale-fact witnesses passes. Provider897function/14target joins, full38/187960words and original C8/all32 pass. Original ABBA22/Decimal/all8/frame results are admitted only if actually completed above; predecessor AL finite results are not transferred. No full science/Mandatory/native-window/installed-release/public-main readiness claim.')
if frame_owner:manifest['frame_terminal']={k:frame_owner.get(k) for k in ('source','passed','wall_seconds','peak_sampled_rss_kib','minimum_sampled_available_kib','stop_reason','child_exit')}
(B/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(destination=str(B.relative_to(R)),terminal_owners=len(owners),new_products=17,all22=all22,all8=all8,original_frame_passed=frame_passed,archives=[dict(file=a['file'],bytes=a['identity']['bytes'],sha256=a['identity']['sha256']) for a in archives])))
