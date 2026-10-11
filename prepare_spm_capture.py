"""Prepare one process-local counter-gate diagnostic on the original control."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
TAG = 'evidence/issue-44-native-executed-cost-6be3132'
assert ROOT == Path('/home/astra/.project/Sirius') and WORK == ROOT/'out/native-spm-cost'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*args): return subprocess.check_output(['git', *args], text=True).strip()
def blob(name): return subprocess.check_output(['git', 'show', TAG+':'+name])
def replace(s, old, new):
    assert s.count(old) == 1, (old[:100], s.count(old))
    return s.replace(old, new)
outputs = {}
def save(name, data):
    assert name not in outputs and not (WORK/name).exists()
    outputs[name] = data.encode() if isinstance(data, str) else data
def saved(name): return outputs[name].decode()

head = read(WORK/'source.json')['source_revision']
assert git('rev-parse', 'HEAD') == head == git('rev-parse', '@{u}')
assert not git('status', '--porcelain')
assert git('rev-parse', 'HEAD^{tree}') == git('rev-parse', 'b81e061^{tree}')
assert read(WORK/'saved-counter-result.json')['counter_data_available'] is False
ancestry = {}
for name in ('owned_console_capture.cs', 'emergency_profile_cleanup.ps1',
             'capture-reference.json', 'tools-extraction.json', 'tools-gui-extraction.json'):
    data = blob(name)
    save(name, data)
    ancestry[name] = {'revision': git('rev-parse', TAG), 'path': name,
                      'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
original = {}
for name in ('run_interval_capture_control.ps1', 'run_native_interval_capture.py',
             'verify_interval_capture_control.py', 'capture-interval-plan.json'):
    data = blob(name)
    original[name] = data.decode()
    ancestry[name] = {'revision': git('rev-parse', TAG), 'path': name,
                      'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

ps = original['run_interval_capture_control.ps1']
ps = replace(ps, '$ErrorActionPreference=\'Stop\'', '''$ErrorActionPreference='Stop'
$counterVariable='RDP_ENABLE_SPM_COUNTERS'
$counterBefore=[Environment]::GetEnvironmentVariable($counterVariable,'Process')
$counterUser=[Environment]::GetEnvironmentVariable($counterVariable,'User')
$counterMachine=[Environment]::GetEnvironmentVariable($counterVariable,'Machine')
if ($null -ne $counterBefore -or $null -ne $counterUser -or $null -ne $counterMachine) {throw 'Unexpected existing counter support-gate environment control'}''')
ps = ps.replace('original-legacy-interval.rgp', 'spm-original-legacy-interval.rgp')
ps = replace(ps, 'No explicit clock or driver experiment configuration.',
             'One diagnostic RDP-child-only support-gate environment override; no user/machine environment, registry, driver installation, scientific source, product or stream change.')
ps = replace(ps, '$record.environment=[ordered]@{}', '''$record.counter_environment=[ordered]@{name=$counterVariable;value='1';before_process=$counterBefore;before_user=$counterUser;before_machine=$counterMachine;set_at_panel_creation=$false;restored_before_panel_resume=$false;absent_before_science=$false;restored_at_closure=$false;scope='Experimental public-source support-gate bypass in the owned RDP child; actual recorded SPM configuration/data must separately establish effectiveness, not hardware qualification.'}
$record.environment=[ordered]@{}''')
ps = replace(ps, "$record.environment[$name]=[Environment]::GetEnvironmentVariable($name) }",
             "$record.environment[$name]=[Environment]::GetEnvironmentVariable($name) }\n$record.environment[$counterVariable]=$counterBefore")
ps = replace(ps, '$consolePid=[SiriusOwnedConsole]::StartSuspended($panelExe,($panelArguments -join \' \'),$Tools)', '''# Only CreateProcessW for the suspended RDP child occurs while the process variable exists.
 # RDS was created before this control; restore before resuming RDP or creating science.
 try {
  [Environment]::SetEnvironmentVariable($counterVariable,'1','Process')
  if ([Environment]::GetEnvironmentVariable($counterVariable,'Process') -ne '1') {throw 'Counter support-gate process setting failed'}
  $consolePid=[SiriusOwnedConsole]::StartSuspended($panelExe,($panelArguments -join ' '),$Tools)
  $record.counter_environment.set_at_panel_creation=$true
 } finally {
  [Environment]::SetEnvironmentVariable($counterVariable,$counterBefore,'Process')
 }
 if ($null -ne [Environment]::GetEnvironmentVariable($counterVariable,'Process')) {throw 'Counter support-gate process environment not restored before panel resume'}
 $record.counter_environment.restored_before_panel_resume=$true''')
ps = replace(ps, '$watch.Restart()\n $process=Start-Process', '''if ($null -ne [Environment]::GetEnvironmentVariable($counterVariable,'Process')) {throw 'Counter support-gate environment leaked before science launch'}
 $record.counter_environment.absent_before_science=$true
 $watch.Restart()
 $process=Start-Process''')
ps = replace(ps, "$cleanupErrors=New-Object 'System.Collections.Generic.List[string]'", '''$cleanupErrors=New-Object 'System.Collections.Generic.List[string]'
 try {
  [Environment]::SetEnvironmentVariable($counterVariable,$counterBefore,'Process')
  $record.counter_environment.restored_at_closure=($null -eq [Environment]::GetEnvironmentVariable($counterVariable,'Process') -and $null -eq [Environment]::GetEnvironmentVariable($counterVariable,'User') -and $null -eq [Environment]::GetEnvironmentVariable($counterVariable,'Machine'))
  if (!$record.counter_environment.restored_at_closure) {throw 'Counter control environment closure mismatch'}
 } catch {$cleanupErrors.Add('counter-environment-restore: '+$_.Exception.Message)}''')
save('run_spm_capture_control.ps1', ps)

reader = original['verify_interval_capture_control.py']
reader = reader.replace('capture-interval-plan.json', 'capture-spm-plan.json')
reader = reader.replace('original-legacy-interval.rgp', 'spm-original-legacy-interval.rgp')
reader = reader.replace('native-capture-interval-console', 'native-spm-capture-console')
reader = replace(reader, "assert not o['accessible_remaining_developer_processes']", '''ce=o['counter_environment']
    assert ce['name']=='RDP_ENABLE_SPM_COUNTERS' and ce['value']=='1'
    for k in ('before_process','before_user','before_machine'):assert ce[k] is None
    for k in ('set_at_panel_creation','restored_before_panel_resume','absent_before_science','restored_at_closure'):assert ce[k] is True
    assert not o['accessible_remaining_developer_processes']''')
save('verify_spm_capture_control.py', reader)

bridge = original['run_native_interval_capture.py']
bridge = bridge.replace('native-capture-interval-console', 'native-spm-capture-console')
bridge = bridge.replace('native-capture-interval-source-review.json', 'native-spm-capture-source-review.json')
bridge = bridge.replace('run_interval_capture_control.ps1', 'run_spm_capture_control.ps1')
bridge = bridge.replace('verify_interval_capture_control', 'verify_spm_capture_control')
bridge = bridge.replace('capture-interval-plan.json', 'capture-spm-plan.json')
bridge = bridge.replace('original-legacy-interval.rgp', 'spm-original-legacy-interval.rgp')
bridge = replace(bridge, "assert read(WORK/'service-info-loopback-bridge.json')['accepted']",
                 "# The unchanged acquired interval ancestry is sealed in the plan; current owned readiness is required again.")
begin = bridge.index('paths=[Path(__file__)')
end = bridge.index('seals={seal_key(p)', begin)
bridge = bridge[:begin]+'''paths=[Path(__file__),WORK/'run_spm_capture_control.ps1',WORK/'emergency_profile_cleanup.ps1',WORK/'source.json',WORK/'tools-extraction.json',WORK/'tools-gui-extraction.json',WORK/'native-spm-capture-source-review.json',WORK/'capture-spm-plan.json',WORK/'capture-reference.json',WORK/'verify_spm_capture_control.py',WORK/'owned_console_capture.cs',WORK/'prepare_spm_capture.py',WORK/'counter-route-review.json',WORK/'saved-counter-result.json']
paths += [Path(p) if Path(p).is_absolute() else ROOT/p for p in plan['files']]
''' + bridge[end:]
bridge = replace(bridge, "assert read(WORK/'console-requested-prefix-result.json')['bounded_readiness_prefix_verified']",
                 "# Same reviewed counted-console adapter; fresh actual width and cell-count checks remain mandatory.")
bridge = replace(bridge, 'seven-dispatch legacy-repeat0 coverage range; per-event sample availability unproved;',
                 'seven-dispatch legacy-repeat0 coverage with one owned-RDP-only SPM support-gate environment control; actual SPM configuration and sample chunks require separate analysis;')
save('run_native_spm_capture.py', bridge)

plan = json.loads(original['capture-interval-plan.json'])
files = {p: v for p, v in plan['files'].items() if p.startswith(('bin/', '/mnt/', 'attestations/'))}
for p, value in files.items(): assert seal(Path(p) if Path(p).is_absolute() else ROOT/p) == value, p
assert len([p for p in files if p.startswith('bin/')]) == 164
for name in ('tools-extraction.json', 'tools-gui-extraction.json'):
    for v in json.loads(saved(name))['files']:
        assert files[v['path']] == {k: v[k] for k in ('bytes', 'sha256')}
prior = ROOT/'attestations/native-vulkan/6be3132/native-executed-cost/manifest.json'
files[str(prior.relative_to(ROOT))] = seal(prior)
for p in WORK.iterdir():
    if p.is_file(): files[str(p.relative_to(ROOT))] = seal(p)
for name, data in outputs.items():
    files[str((WORK/name).relative_to(ROOT))] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
dll = ROOT/plan['tools']/'RdpCaptureApi.dll'
assert b'RDP_ENABLE_SPM_COUNTERS' in dll.read_bytes()
plan.update(files=files, ancestral_inputs=ancestry,
    diagnostic_control={'name':'RDP_ENABLE_SPM_COUNTERS','value':'1',
        'scope':'Only owner process while creating suspended owned RDP child; parent restored before RDP resume and scientific child creation; never user/machine environment.',
        'reference_revision':'d63073a854c8ccc4b10edbf27559bb698c5685b3',
        'recorded_asic_family':148,'recorded_driver_base':'26.10.41.01',
        'actual_capture_api_literal_present':True,'capture_api':seal(dll),
        'limits':'Public-source route and packaged literal do not attest hardware support or binary internal eligibility. Actual trace must show enabled configuration, nonempty requested counters and nonempty SPM sample chunks.'},
    scope='One changed counter-gate control on the unchanged previously acquired dispatch:9:7 original numerical test, 90s/4GiB science and existing finite tool guards. Fresh complete buffers, seven modules, 51 markers, actual provider and trace required. No additional dispatches, new shader instrumentation, clock/driver/global-environment change, performance benefit, rendered frame or qualification claim.')
save('capture-spm-plan.json', json.dumps(plan, indent=2)+'\n')
for name, data in outputs.items():
    if name.endswith('.py'): ast.parse(data, filename=name)
for name, data in outputs.items(): (WORK/name).write_bytes(data)
print(json.dumps({'prepared_files': len(outputs), 'fresh_plan_input_versions': len(files),
                  'ancestral_revision': git('rev-parse', TAG), 'native_and_tools':164}))
