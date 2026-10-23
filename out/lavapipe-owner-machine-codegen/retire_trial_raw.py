"""Remove only pushed archive-verified completed raw output and unused staging."""
from pathlib import Path
import hashlib, json, os, subprocess, zipfile
R=Path.cwd(); W=Path(__file__).resolve().parent
D=R/'attestations/qualification/33623035/lavapipe-owner-machine-codegen-finite-outcome'
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
assert not (D/'cleanup.json').exists()
remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/development/lavapipe-owner-machine-codegen-2026-10-23'],text=True).split()[0]
assert subprocess.run(['git','merge-base','--is-ancestor','ec993a02',remote]).returncode==0
manifest=doc(D/'MANIFEST.json');closure=doc(D/'ownership-closure.json')
assert closure['all_owned_runs_terminal_and_reaped']
for owner in closure['owners'].values():
 for birth in owner['observed_births']:
  p=Path('/proc')/str(birth['pid'])/'stat'
  if p.exists():assert p.read_text().rsplit(')',1)[1].split()[19]!=birth['saved_start_ticks']
for name,expected in manifest['archives'].items():assert pin(D/name)==expected
tracked={R/os.fsdecode(p) for p in subprocess.check_output(['git','ls-files','-z']).split(b'\0') if p}
targets={}
with zipfile.ZipFile(D/'diagnostics-and-trial.zip') as z:
 for key,expected in manifest['diagnostics'].items():
  p=R/key
  if p in tracked or not p.is_file():continue
  # Keep current configured-build proofs and concise owners/code findings.
  if not (p.is_relative_to(W/'ordinary') or p.is_relative_to(W/'runtime') or p.is_relative_to(W/'readers')):continue
  if p.name in {'owner.json','capture-verification.json','ordinary-code-review.json','benefit.json','decimal-agreement.json','control-8-acceptance.json','sole-ABBA-inputs.json'}:continue
  assert pin(p)=={k:expected[k] for k in ['bytes','sha256']}
  b=z.read(key);assert len(b)==expected['bytes'] and hashlib.sha256(b).hexdigest()==expected['sha256']
  targets[p]=expected['bytes']
# This rejected staged copy has no runtime consumer. The configured build's
# actual target, exact baseline and mixed-host products stay available.
for key,expected in manifest['actual_product_paths'].items():
 if not key.startswith('bin/linux-gcc/lavapipe-owner-machine-codegen/provider/'):continue
 p=R/key;assert pin(p)=={k:expected[k] for k in ['bytes','sha256']}
 blob=manifest['content_addressed_blobs'][expected['sha256']]
 archive=R/blob['archive'] if blob.get('reused') else D/blob['archive']
 with zipfile.ZipFile(archive) as z:
  b=z.read(blob['member']);assert len(b)==expected['bytes'] and hashlib.sha256(b).hexdigest()==expected['sha256']
 targets[p]=expected['bytes']
publication=doc(D/'publication.json');receipt=D/'publication-receipts.zip'
assert pin(receipt)==publication['receipt_archive']
with zipfile.ZipFile(receipt) as z:
 for key,expected in publication['receipt_members'].items():
  p=R/key
  if p in tracked or not p.is_file():continue
  assert pin(p)==expected and z.read(key)==p.read_bytes()
  targets[p]=expected['bytes']
for p in targets:p.unlink()
for folder in [W/'ordinary',W/'runtime',W/'readers',R/'bin/linux-gcc/lavapipe-owner-machine-codegen/provider']:
 if folder.exists():
  for p in sorted(folder.rglob('*'),key=lambda p:len(p.parts),reverse=True):
   if p.is_dir() and not any(p.iterdir()):p.rmdir()
  if not any(folder.iterdir()):folder.rmdir()
result={'preserved_and_pushed_before_deletion':'ec993a02','all10_owned_runs_terminal_and_saved_births_absent':True,
 'removed_files':len(targets),'removed_bytes':sum(targets.values()),
 'retained':['All tracked source and pushed finite bundles','Concise owner, code and original-gate findings; completed owner folders continue to refuse unchanged reruns',
 'Configured Mesa source/build and sealed build proofs for ongoing latency development','Original baseline, immutable modules, complete-input products and mixed-host products'],
 'scope':'Only authentic archive-verified completed raw captures/runtime/reader output, duplicate publication copies and unused rejected staging. No source/build-context deletion, acceptance retry or release-admission inference.'}
(D/'cleanup.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
