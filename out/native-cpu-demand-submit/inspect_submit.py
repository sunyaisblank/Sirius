"""Check the changed queue path in the actual provider without running a workload."""
from pathlib import Path
import hashlib,json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
b=json.loads((W/'provider/provider-build-result.json').read_text());p=R/b['library']
assert b['passed'] and p.stat().st_size==b['library_identity']['bytes']
assert hashlib.file_digest(p.open('rb'),'sha256').hexdigest()==b['library_identity']['sha256']
code={name:subprocess.check_output(['objdump','-d','--disassemble='+name,str(p)],text=True)
      for name in ('lvp_CreateDevice','vk_queue_submit')}
for name,text in code.items():assert '<'+name+'>:' in text
assert '<vk_queue_init>' in code['lvp_CreateDevice']
assert '<vk_queue_enable_submit_thread' not in code['lvp_CreateDevice']
targets=sorted(set(re.findall(r'call\s+[^<]+<([^>]+)>',code['vk_queue_submit'])))
assert {'vk_sync_wait_many','thrd_create','vk_queue_submit_final','vk_queue_submit_move_binary_waits_to_temps','vk_sync_reset'}<=set(targets)
current=json.loads((W/'provider/source-contract.json').read_text())
old=json.loads((R/'out/native-cpu-masked-arithmetic/provider/source-contract.json').read_text())
assert current['native_objects']==old['native_objects']
result={'passed':True,'provider':b['library'],'provider_identity':b['library_identity'],
 'device_init_forced_submit_thread_calls':0,'runtime_submit_call_targets':targets,
 'all8_native_objects_unchanged':True,
 'scope':'Actual changed-provider code check: device creation calls queue init without forced thread startup; runtime submit retains pending waits, thread creation, final submission and binary reset paths. Exact source review supplies dependency/ordering reasoning. This is not a machine-checked control-flow proof, executed unresolved-wait test, performance or scientific/release acceptance.'}
with (W/'submit-code-inspection.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'passed':True,'native_objects_unchanged':8,'forced_submit_thread_calls':0,'runtime_thread_fallback_retained':True}))
