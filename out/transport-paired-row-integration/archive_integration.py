"""Preserve this completed finite native trial before restoration or cleanup."""
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
source=json.loads((WORK/'source.json').read_text());head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def identity(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
controls=read(WORK/'native-controls-analysis.json');assert controls['pass'] and controls['total_cases']==10
comparison=read(WORK/'native-retention-comparison.json');assert len(comparison['checks'])==46
assert comparison['candidate_revision']==head and comparison['numerical_stream_provider_owner_prerequisites_pass']
assert read(WORK/'native-matched-sequence.json')['completed']
preservation=read(WORK/'source-git-preservation.json');assert preservation['pushed_verified'] and preservation['base']==head
live_tools={str(p.relative_to(ROOT)) for p in WORK.rglob('*') if p.is_file() and p.suffix in ('.py','.ps1') and '__pycache__' not in p.parts and 'temp' not in p.relative_to(WORK).parts}
assert live_tools=={r['path'] for r in preservation['sources']}
for r in preservation['sources']:
 assert identity(ROOT/r['path'])=={k:r[k] for k in ('bytes','sha256')}
 data=subprocess.check_output(['git','show',preservation['commit']+':'+r['path']])
 assert {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}=={k:r[k] for k in ('bytes','sha256')}
DEST=ROOT/'attestations/native-vulkan'/head[:7]/'paired-transport-integrated-trial'
assert not DEST.exists()
bindings=read(WORK/'execution-bindings.json')['files']
for relative,record in bindings.items():
 p=Path(relative);p=p if p.is_absolute() else ROOT/p
 assert identity(p)==record,relative
software={}
for r in read(WORK/'controls-complete.json')['results']:
 owner=read(WORK/'controls'/r['name']/'owner.json')
 for path,record in owner['whole_input_seals'].items():
  assert path not in software or software[path]==record;software[path]=record
for path,record in software.items():assert identity(Path(path))==record,path
zips=set();retained_exports={}
for label in ('baseline','candidate'):
 for platform in ('windows','macos'):
  receipt=read(WORK/(label+'-'+platform+'-export-download.json'))
  p=WORK/(label+'-'+platform+'-export.zip');zips.add(p)
  assert identity(p)=={'bytes':receipt['download_bytes'],'sha256':receipt['sha256']}
  assert receipt['artifact']['digest']=='sha256:'+receipt['sha256']
  bundle=ROOT/'attestations/native-build'/receipt['source_revision'][:7]/(platform+'-build')
  assert {str(p.relative_to(bundle)) for p in bundle.rglob('*') if p.is_file()}=={r['path'] for r in receipt['extracted_files']}
  for r in receipt['extracted_files']:
   p=bundle/r['path'];record={k:r[k] for k in ('bytes','sha256')};assert identity(p)==record
   retained_exports[str(p.relative_to(ROOT))]=record
selected={};external=dict(retained_exports)
def select(p,target):
 assert p.is_file() and not p.is_symlink()
 assert target not in selected or selected[target]==p;selected[target]=p
for p in sorted(WORK.rglob('*')):
 assert not p.is_symlink()
 if p.is_file() and '__pycache__' not in p.parts and p not in zips and p.name not in ('source-git.index','archive-protection.json'):
  select(p,'diagnostic/'+str(p.relative_to(WORK)))
for producer in read(WORK/'execution-bindings.json')['producers'].values():
 stage=ROOT/producer['stage']
 for p in sorted(stage.rglob('*')):
  assert not p.is_symlink()
  if p.is_file():select(p,'native-stage-inputs/'+str(p.relative_to(ROOT)))
for relative in set(bindings)|set(software):
 p=Path(relative);p=p if p.is_absolute() else ROOT/p
 if not p.is_relative_to(ROOT):continue
 if p.is_relative_to(WORK) or any(p==q for q in selected.values()):continue
 if p.is_relative_to(ROOT/'attestations'):
  external[str(p.relative_to(ROOT))]=identity(p);continue
 select(p,'source-build-inputs/'+str(p.relative_to(ROOT)))
payloads=[]
for target,p in sorted(selected.items()):
 original=identity(p);q=DEST/target;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
 assert identity(q)==original
 payloads.append({'path':target,'origin':str(p.relative_to(ROOT)),**original})
assert {str(p.relative_to(DEST)) for p in DEST.rglob('*') if p.is_file()}=={r['path'] for r in payloads}
for r in payloads:
 expected={k:r[k] for k in ('bytes','sha256')}
 assert identity(DEST/r['path'])==expected and identity(ROOT/r['origin'])==expected
for relative,record in bindings.items():
 p=Path(relative);p=p if p.is_absolute() else ROOT/p
 assert identity(p)==record,relative
for path,record in software.items():assert identity(Path(path))==record,path
for relative,record in external.items():assert identity(ROOT/relative)==record,relative
for r in preservation['sources']:assert identity(ROOT/r['path'])=={k:r[k] for k in ('bytes','sha256')}
manifest={'candidate_revision':head,'baseline_revision':source['baseline_revision'],
 'retention_gate_pass':comparison['retention_gate_pass'],'failed_checks':len(comparison['failed_checks']),
 'payloads':payloads,'payload_count':len(payloads),'payload_bytes':sum(p['bytes'] for p in payloads),
 'retained_external_payloads':[{'path':p,**r} for p,r in sorted(external.items())],
 'scope':'Finite isolated/software/native integration and one unchanged frozen comparison; all original failures retained. No original frame, full scientific/provider/release or main-admission claim.'}
(DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
receipt={'archive':str(DEST.relative_to(ROOT)),'manifest':identity(DEST/'manifest.json'),
 'payload_count':manifest['payload_count'],'payload_bytes':manifest['payload_bytes'],
 'retained_external_payloads':len(external),'retention_gate_pass':comparison['retention_gate_pass']}
(WORK/'archive-protection.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
