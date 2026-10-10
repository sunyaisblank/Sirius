"""Preserve verified restoration and trial closure before finished scratch removal."""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent;OLD=ROOT/'out/transport-paired-row-integration'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def identity(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
source=read(WORK/'source.json');head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
assert subprocess.check_output(['git','rev-parse','HEAD^{tree}'],text=True).strip()==source['source_tree']
assert source['entire_tracked_tree_exact_baseline']
assert read(WORK/'build-and-arrays.json')['restored_all_40_arrays_exact']
assert read(WORK/'linux-readonly-payloads.json')['pass_'] and read(WORK/'build-terminal-groups.json')['pass']
controls=read(WORK/'controls-complete.json');assert controls['returncode']==0 and len(controls['results'])==8
assert read(WORK/'independent-restoration-evidence-review.json')['pass']
final=read(OLD/'final-native-cleanup-bridge.json');assert final['accepted'] and final['source_and_input_seals_exact_at_end']
ci=read(WORK/'restoration-ci.json');assert ci['status']=='completed' and ci['conclusion']=='success' and ci['headSha']==head
assert len([j for j in ci['jobs'] if j['conclusion']=='success'])==3 and len([j for j in ci['jobs'] if j['conclusion']=='skipped'])==4
protection=read(OLD/'archive-protection.json');primary=ROOT/protection['archive'];assert identity(primary/'manifest.json')==protection['manifest']
prior=read(primary/'manifest.json');prior_origins={r['origin']:{k:r[k] for k in ('bytes','sha256')} for r in prior['payloads']}
preserved=read(WORK/'source-git-preservation.json');assert preserved['pushed_verified'] and preserved['base']==head
tools={str(p.relative_to(ROOT)) for p in WORK.rglob('*') if p.is_file() and p.suffix in ('.py','.ps1') and '__pycache__' not in p.parts and 'temp' not in p.relative_to(WORK).parts}
assert tools=={r['path'] for r in preserved['sources']}
for r in preserved['sources']:
 data=subprocess.check_output(['git','show',preserved['commit']+':'+r['path']]);expected={k:r[k] for k in ('bytes','sha256')}
 assert identity(ROOT/r['path'])==expected and {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}==expected
seals={}
for result in controls['results']:
 owner=read(WORK/'controls'/result['name']/'owner.json')
 assert owner['source_and_inputs_unchanged'] and owner['owned_birth_absent'] and owner['owned_group_absent']
 for path,r in owner['whole_input_seals'].items():assert path not in seals or seals[path]==r;seals[path]=r
for path,r in seals.items():assert identity(Path(path))==r,path
DEST=ROOT/'attestations/native-vulkan'/head[:7]/'paired-transport-restoration';assert not DEST.exists()
selected={}
def select(p,target):assert p.is_file() and not p.is_symlink();assert target not in selected or selected[target]==p;selected[target]=p
for p in WORK.rglob('*'):
 assert not p.is_symlink()
 if p.is_file() and '__pycache__' not in p.parts and p.name not in ('archive-protection.json','source-git.index'):select(p,'diagnostic/'+str(p.relative_to(WORK)))
for p in OLD.rglob('*'):
 assert not p.is_symlink()
 if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.zip':
  if prior_origins.get(str(p.relative_to(ROOT)))!=identity(p):select(p,'trial-closure/'+str(p.relative_to(OLD)))
for path in seals:
 p=Path(path)
 if p.is_relative_to(ROOT) and not p.is_relative_to(WORK):select(p,'source-build-inputs/'+str(p.relative_to(ROOT)))
payloads=[]
for target,p in sorted(selected.items()):
 expected=identity(p);q=DEST/target;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);assert identity(q)==expected
 payloads.append({'path':target,'origin':str(p.relative_to(ROOT)),**expected})
assert {str(p.relative_to(DEST)) for p in DEST.rglob('*') if p.is_file()}=={r['path'] for r in payloads}
for r in payloads:
 expected={k:r[k] for k in ('bytes','sha256')};assert identity(DEST/r['path'])==expected and identity(ROOT/r['origin'])==expected
for path,r in seals.items():assert identity(Path(path))==r,path
assert tools=={str(p.relative_to(ROOT)) for p in WORK.rglob('*') if p.is_file() and p.suffix in ('.py','.ps1') and '__pycache__' not in p.parts and 'temp' not in p.relative_to(WORK).parts}
for r in preserved['sources']:assert identity(ROOT/r['path'])=={k:r[k] for k in ('bytes','sha256')}
manifest={'source_revision':head,'source_tree':source['source_tree'],'exact_baseline_restoration':True,'original_controls_passed':8,'payloads':payloads,'payload_count':len(payloads),'payload_bytes':sum(r['bytes'] for r in payloads),'primary_trial_archive':protection,'scope':'Strict exact-source restoration/build/all40arrays/160ELFjoins/eightoriginalcontrols/threeplatformnonrenderCI/fresh45birthterminalclosure after rejected one46checkcomparison. No originalframe/fullscientific/provider/release/mainadmission claim.'}
(DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'archive':str(DEST.relative_to(ROOT)),'manifest':identity(DEST/'manifest.json'),'payload_count':manifest['payload_count'],'payload_bytes':manifest['payload_bytes']}
(WORK/'archive-protection.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
