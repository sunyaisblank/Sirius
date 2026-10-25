"""Seal one execution of the explicitly admitted, previously completed link."""
from pathlib import Path
import hashlib,json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'fixture';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
assert not subprocess.check_output(['git','status','--porcelain']).strip()
p=J(F/'prepared.json');a=J(F/'build-admission.json');assert a['passed'] and a['original_owner_passed'] is False
assert pin(F/'build/owner.json')==a['original_owner'];assert not (F/'run').exists()
inputs=set(p['inputs'])|set(a['runtime_libraries'])|{p['executable']}
for name in ['fixture/prepared.json','fixture/build-result.json','fixture/build/owner.json','fixture/build-admission.json','fixture/linked-code-admission.json','admit_fixture_build.py','prepare_fixture_run.py','run_fixture.py']:
 inputs.add(str((W/name).relative_to(R)))
for name in inputs:assert (R/name).is_file(),name
actions=J(F/'actions.json');assert 'run' not in actions
actions['run']=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
 limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=['/usr/bin/python3','-B',str(W/'run_fixture.py')],frozen_inputs=sorted(inputs),
 scope=p['scope']+' Original build owner remains failed at its reader; an explicit zero-recompile/link admission binds the existing executable and all three objects/GOT targets. One first numerical workload only.')
(F/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps(dict(prepared='run',repeated_builds=0)))
