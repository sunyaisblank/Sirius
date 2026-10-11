"""Remove only the sealed superseded diagnostic cache after fresh scoped checks."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, os, shutil, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
STAGE = ROOT / 'bin/windows-msvc/native-8f1118d-shader-clock-stream'
ARCHIVE = ROOT / 'attestations/native-vulkan/8f1118d/native-shader-clock-execution'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    data = p.read_bytes()
    return {'bytes':len(data), 'sha256':hashlib.sha256(data).hexdigest()}
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip() == '24a8e145d1da5dca75a4fa26ec8e639020eecf37'
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
assert read(WORK / 'results.json')['accepted'] and read(WORK / 'actual-results-review.json')['accepted']
assert read(WORK / 'obsolete-stage-cleanup/retirement-run.json')['accepted']
census = read(WORK / 'obsolete-stage-cleanup/native-retirement.json')
assert census['pass'] and not census['matches'] and not census['unused_stage_matches']
assert census['unused_stage'].lower() == subprocess.check_output(['wslpath','-w',str(STAGE)],text=True).strip().lower()
assert seal(ARCHIVE / 'manifest.json')['sha256'] == 'e9bc1c10ea08c866f124cca6420c6896c82e1d23765bc747c05a7c82d74ca208'
plan = read(ARCHIVE / 'files/execution-plan.json')
prefix = str(STAGE.relative_to(ROOT)) + '/'
expected = {k:v for k,v in plan['input_seals'].items() if k.startswith(prefix)}
assert len(expected) == 83 and STAGE.is_dir() and not STAGE.is_symlink()
entries = list(STAGE.rglob('*'))
assert all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in entries)
actual = {str(p.relative_to(ROOT)):seal(p) for p in entries if p.is_file()}
assert actual == expected
payloads = {(v['bytes'],v['sha256']):str(p.relative_to(ROOT)) for p in (ARCHIVE / 'files/artifact').rglob('*') if p.is_file() for v in [seal(p)]}
assert all((v['bytes'],v['sha256']) in payloads for v in expected.values())
# Check accessible Linux command lines separately; no hidden/scheduled claim.
matches = []
for d in Path('/proc').iterdir():
    if not d.name.isdigit() or int(d.name) == os.getpid(): continue
    try: command = (d/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
    except (FileNotFoundError,ProcessLookupError,PermissionError): continue
    if str(STAGE) in command or prefix in command: matches.append({'pid':int(d.name),'command':command})
assert not matches, matches
receipt = {'removed_stage':prefix.rstrip('/'), 'files':83, 'bytes':sum(v['bytes'] for v in expected.values()),
           'recovery_payloads':len({(v['bytes'],v['sha256']) for v in expected.values()}),
           'fresh_windows_census':seal(WORK/'obsolete-stage-cleanup/native-retirement.json'),
           'linux_accessible_matches':matches, 'scope':'Exact sealed obsolete cache, whole-byte recovery in closed archive; fresh accessible scoped checks. No hidden/scheduled-consumer absence claim.'}
shutil.rmtree(STAGE)
assert not STAGE.exists()
receipt['removed'] = True
(WORK/'obsolete-stage-removal.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
