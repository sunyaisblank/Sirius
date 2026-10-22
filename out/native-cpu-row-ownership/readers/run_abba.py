"""Exactly one ordered trial; each existing owned slot refuses an existing folder."""
from pathlib import Path
import hashlib,json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;RUN=R/'out/native-cpu-row-ownership/runtime'
acceptance=W/'control-8-acceptance.json';assert json.loads(acceptance.read_text())['passed']
assert not any((RUN/label).exists() for label in ['a1','b1','b2','a2'])
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
inputs={str(p.relative_to(R)):pin(p) for p in [*W.glob('*.py'),acceptance,W/'expected-readback-layout.json',RUN/'actions.json',RUN/'run_case.py',RUN/'run_owned.py',RUN/'owned_guard.py']}
(W/'sole-ABBA-inputs.json').write_text(json.dumps({'order':['a1','b1','b2','a2'],'inputs':inputs,'extra_warmups_or_resampling':0},indent=2)+'\n')
for label in ['a1','b1','b2','a2']:
 subprocess.run(['/usr/bin/python3','-B',str(RUN/'run_owned.py'),label],check=True)
 assert inputs=={key:pin(R/key) for key in inputs}
subprocess.run(['/usr/bin/python3','-B',str(W/'read_results.py'),'benefit'],check=True)
