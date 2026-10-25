"""Preserve the exact integer identity-boundary component without invented gates."""
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
f=J(W/'fixture/result.json');assert f['passed'] and f['calls']==555660 and f['changed_operator_calls']==333396 and f['compared_words']==66679200 and f['tag_words']==13335840
assert J(W/'identity-code-inspection.json')['passed'] and J(C/'control-8-acceptance.json')['passed'] and J(D/'result.json')['passed']
assert not (C/'benefit').exists() and not (C/'controls').exists() and not (C/'cold-frame').exists()
owners={};boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
for root in (W,C,D):
 for path in sorted(root.rglob('owner.json')):
  o=J(path);assert o['status']=='terminal' and o['passed'] and o['child_exit']==0 and o['source_unchanged'] and o['whole_inputs_unchanged']
  assert not o['cleanup_errors'] and not o['remaining_owned_processes'] and all(x['absent'] for x in o['observed_births'])
  for row in o['observed_births']:
   if boot==o['boot_id']:
    try:t=(Path('/proc')/str(row['pid'])/'stat').read_text().rsplit(')',1)[1].split()[19]
    except (FileNotFoundError,ProcessLookupError):t=None
    assert t!=str(row['saved_start_ticks'])
  owners[str(path.relative_to(R))]={k:o.get(k) for k in ('source','status','passed','child_exit','wall_seconds','peak_sampled_rss_kib','stop_reason')}
assert len(owners)==8
source=J(C/'control-8/owner.json')['source']['revision'];B=R/'attestations/qualification'/source[:8]/'native-cpu-local-identities-finite-compatibility';assert not B.exists();B.mkdir(parents=True)
compiled=J(W/'ispc-dense-compile-result.json');provider=J(W/'provider/provider-build-result.json');fixture=J(W/'fixture/build-result.json');diagnostic=J(D/'link-result.json')
products=[R/x['object']['path'] for x in compiled['objects']]+[R/x['product']['path'] for x in fixture['commands']]+[R/provider[k] for k in ('library','ICD')]+[R/diagnostic['executable']['path']]
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
archives=[archive('source-build-and-compatibility.zip',sources),archive('actual-products.zip',{str(p.relative_to(R)):p for p in products})]
reference=R/'attestations/qualification/d4f1d31e/native-cpu-value-fact-stores-frame-timeout/MANIFEST.json';assert subprocess.check_output(['git','show','@{upstream}:'+str(reference.relative_to(R))])==reference.read_bytes()
manifest=dict(source=source,accepted=False,finite_component_passed=True,original_C8_passed=True,full38_word_comparison_passed=True,original_ABBA22_attempted=False,original_all8_controls_completed=False,original_frame_attempted=False,full_scientific_and_release_qualification=False,new_products=9,owners=owners,archives=archives,production_raw_compilations=3,production_CPP_compilations=0,diagnostic_ISPC_compilations=1,diagnostic_CPP_compilations=0,diagnostic_symbol_renames=1,diagnostic_links=2,reused_provenance={str(reference.relative_to(R)):identity(reference)},
 scope='AN moves the original AM complete-fact Add/Subtract/Multiply identity/refusal stores into a common forced-inline integer-only ISPC body before the same C++ ABI calls. All active lanes must be complete and choose an original identity/refusal; otherwise return before any payload store and enter unchanged AM math. Copies retain noncanonical nonzero valid/signed-zero rules; negation toggles only original three sign bits; zero/refusal constants and metadata match AM. Snapshot selected operands before aliasing stores. Same context/bank/tag arrays, all writer bodies,9380canonical tuples, general arithmetic/Divide/Sqrt refinements, guards, FP/ISA/count/cache/governors/scheduler/fields/RHS/trials/quartics/rollback remain. No new persistent allocation or floating arithmetic. Actual common-body fixture uses unchanged AM math and original driver:333396changed operators among555660calls,66679200bank/13335840metadata words and1666980active-lane results pass, including fresh/unknown/legacy metadata. Driver allocated bytes and all inverse-normalized symbols/relocations match after only three undefined-name changes; linked fixture/fallback targets join. All nine actual raw wrappers contain integer-only early-return branches and original fallback. Provider897functions/14targets/refinement joins, full38/187960words and C8/all32/51markers/54joins pass. No runtime hit count or gain is claimed. This is a small component awaiting materially justified whole-path adoption work; no new ABBA/all8/frame cycle was planned or run, and predecessor AM results are not transferred. Public source remains0109ab46; no full science/Mandatory/native-window/installed-release/exact-main readiness.')
(B/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(destination=str(B.relative_to(R)),terminal_owners=len(owners),new_products=9,archived_members=sum(len(a['members']) for a in archives),archives=[dict(file=a['file'],bytes=a['identity']['bytes'],sha256=a['identity']['sha256']) for a in archives])))
