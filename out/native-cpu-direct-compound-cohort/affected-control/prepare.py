"""Run original connected-interval assertions on the reused host and changed provider once."""
from pathlib import Path
import hashlib,json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;C=W.parent
D=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
assert not subprocess.check_output(['git','status','--porcelain']).strip()
assert not (W/'actions.json').exists() and not (W/'control-5').exists()
assert D(C/'control-8-acceptance.json')['passed']
c8=D(C/'control-8/owner.json');assert c8['passed'] and c8['source_unchanged'] and c8['whole_inputs_unchanged']
assert not c8['remaining_owned_processes'] and not c8['cleanup_errors']
h=R/'out/native-cpu-owned-lower-increments';host=D(h/'host-build-result.json')
assert host['passed'] and host['new_cpp_compilations']==1
provider=D(R/'out/native-cpu-direct-compound/provider/provider-build-result.json');assert provider['passed']
exe=R/host['executable'];assert pin(exe)==host['executable_identity']
original=D(R/'out/native-original-input-ownership/runtime/actions.json')['control-5']
assert original['case']=='RetainedDopriTest.ConnectedIntervalsPreserveTrialsAndIndependentBudgets'
assert original['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
inputs=set(D(C/'actions.json')['control-8']['frozen_inputs'])
inputs.update(host['source_inputs']);inputs.update(host['runtime_libraries'])
inputs.update(str(p.relative_to(R)) for p in W.iterdir() if p.is_file())
for p in (exe,h/'host-build-result.json',h/'build/host-build/owner.json',C/'control-8-acceptance.json',C/'control-8/owner.json',C/'control-8/gtest.xml',R/'out/native-original-input-ownership/runtime/actions.json'):
 inputs.add(str(p.relative_to(R)))
assert all((R/name).is_file() for name in inputs)
action={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'case':original['case'],
 'kind':original['kind'],'expected_exit':0,'limits':original['limits'],
 'argv':[str(exe),'--gtest_filter='+original['case'],'--gtest_color=no'],
 'frozen_inputs':sorted(inputs),'software_provider':True,
 'runtime_provider':str(R/provider['library']),'runtime_icd':str(R/provider['ICD']),
 'scope':'One fresh original control5 for reused owned-lower-increment host and AI direct compound opcodes, after current full38/C8. Unmodified original BothProductModes assertions, inputs, independent trials/budgets, endpoint/increment/quarter packets and refusal checks. Original300s/4GiB/2GiB/cache-off. Z broader controls/ABBA remain predecessor evidence, not fresh current qualification. No frame/science/release admission.'}
(W/'actions.json').write_text(json.dumps({'control-5':action},indent=2)+'\n')
print(json.dumps({'prepared':'original-control-5','executed':False}))
