"""Prepare source-only final scoped retirement tools from the reviewed trial controller."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json
ROOT = Path.cwd(); WORK = Path(__file__).resolve().parent
assert WORK == ROOT/'out/native-division-first-correction-screen'
OLD = ROOT/'attestations/native-vulkan/99a4a6d/native-bound-selection-trial/files'
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def put(p, value):
    assert not p.exists(), p
    p.write_text(value)
original = (OLD/'run_trial_retirement.py').read_text()
controller = original.replace('native-99a4a6d-bound-select-trial', 'native-ea6bbb0-division-first-correction-trial').replace('verify_trial_retirement.ps1', 'verify_final_retirement.ps1')
before = "actions=read(WORK/'native-cases.json')['actions']\ninput_paths=(Path(__file__),HERE/'verify_final_retirement.ps1',WORK/'native-cases.json')+tuple(WORK/action/name for action in actions for name in ('run/owner.json','owner-temp/bootstrap.json'))"
after = """review=read(WORK/'closing-retirement-source-review.json')
assert review['accepted']
assert review['controller_seal']=={'bytes':Path(__file__).stat().st_size,'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
assert review['powershell_seal']=={'bytes':(HERE/'verify_final_retirement.ps1').stat().st_size,'sha256':hashlib.sha256((HERE/'verify_final_retirement.ps1').read_bytes()).hexdigest()}
targets=read(WORK/'closing-births.json')
assert targets['accepted'] and len(targets['windows_targets'])==51 and len(targets['linux_targets'])==62
input_paths=(Path(__file__),HERE/'verify_final_retirement.ps1',WORK/'closing-retirement-source-review.json',WORK/'closing-births.json')+tuple(ROOT/p for p in targets['input_seals'])
for p,expected in targets['input_seals'].items():
 assert {'bytes':(ROOT/p).stat().st_size,'sha256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest()}==expected,p"""
assert controller.count(before) == 1
controller = controller.replace(before, after)
controller = controller.replace("WORK/'retirement-bootstrap.json'", "WORK/'final-retirement-bootstrap.json'").replace("WORK/'native-retirement.json'", "WORK/'final-native-retirement.json'").replace("WORK/'retirement.stdout'", "WORK/'final-retirement.stdout'").replace("WORK/'retirement.stderr'", "WORK/'final-retirement.stderr'").replace("WORK/'retirement-run.json'", "WORK/'final-retirement-run.json'")
anchor = "\nargs=[ps,'-NoProfile'"
assert controller.count(anchor) == 1
observe = """
def observe_saved_linux_births():
 observations=[]
 for target in targets['linux_targets']:
  absent=observed_absent(int(target['pid']),str(target['start_ticks']))
  observations.append({**target,'recorded_birth_absent':absent})
  assert absent,('Recorded Linux birth remains live',target)
 return observations
"""
controller = controller.replace(anchor, observe+anchor)
anchor = "try:\n with (WORK/'final-retirement.stdout')"
assert controller.count(anchor) == 1
controller = controller.replace(anchor, "try:\n record['saved_linux_births_before']=observe_saved_linux_births()\n with (WORK/'final-retirement.stdout')")
anchor = "try:\n record['inputs_unchanged']=seal_inputs()==record['input_seals']"
assert controller.count(anchor) == 1
controller = controller.replace(anchor, "try:\n record['saved_linux_births_after']=observe_saved_linux_births()\n record['inputs_unchanged']=seal_inputs()==record['input_seals']")
controller = controller.replace("Fresh exact native control births + accessible scoped consumer census only.", "Fresh exact51 Windows and62 Linux saved births plus accessible scoped WORK/native-stage consumers; no global/hidden/group absence claim.")
put(WORK/'run_final_retirement.py', controller)
original_ps = (OLD/'verify_trial_retirement.ps1').read_text()
start = original_ps.index(' $targets=@()\n')
end = original_ps.index(' foreach($target in $targets){', start)
replacement = """ $births=Get-Content -Raw -LiteralPath (Join-Path $Stage 'closing-births.json')|ConvertFrom-Json
 if(!$births.accepted){throw 'Unaccepted closing birth inventory'}
 $targets=@($births.windows_targets)
 if($targets.Count -ne 51){throw 'Expected34scientificand17auditbirths'}
 if(@($targets|Where-Object {$_.category -eq 'scientific'}).Count -ne 34){throw 'Expected34scientificbirths'}
 if(@($targets|Where-Object {$_.category -eq 'audit'}).Count -ne 17){throw 'Expected17auditbirths'}
 if(@($targets|ForEach-Object {([string]$_.pid)+'|'+$_.birth}|Select-Object -Unique).Count -ne 51){throw 'Duplicatebirth'}
"""
ps = original_ps[:start]+replacement+original_ps[end:]
ps = ps.replace("kind='finite-native-clock-control-retirement'", "kind='finite-first-correction-trial-retirement'")
put(WORK/'verify_final_retirement.ps1', ps)
put(WORK/'closing-retirement-preparation.json', json.dumps({'accepted': True,
    'input_seals': {str(p.relative_to(ROOT)): seal(p) for p in (Path(__file__),OLD/'run_trial_retirement.py',OLD/'verify_trial_retirement.ps1')},
    'output_seals': {n: seal(WORK/n) for n in ('run_final_retirement.py','verify_final_retirement.ps1')},
    'scope': 'Source preparation only: original40s observer and scoped cleanup retained. Final inventory requires all16 original runtime actions plus private query and all17 audit observers.'}, indent=2)+'\n')
print(json.dumps({'accepted': True, 'outputs': {n: seal(WORK/n) for n in ('run_final_retirement.py','verify_final_retirement.ps1')}}))
