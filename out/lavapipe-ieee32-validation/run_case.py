from pathlib import Path
import hashlib,importlib.util,json,os,resource,struct,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-ieee32-staged';label,folder=sys.argv[1:];folder=Path(folder)
action=json.loads((W/'actions.json').read_text())[label];arm=action['arm'];poison=action['poison']
module=R/'bin/linux-gcc/lavapipe-native-fp32/native.spv' if arm=='baseline' else B/'native-preserve.spv'
executable=R/('bin/linux-gcc/lavapipe-fpstate-contract/gpu_probe' if poison else 'bin/linux-gcc/lavapipe-native-fp32/gpu_probe')
driver=B/arm/'libvulkan_lvp.so';oracle=R/'out/lavapipe-native-fp32'
resource.setrlimit(resource.RLIMIT_CORE,(0,0));_,hard=resource.getrlimit(resource.RLIMIT_FSIZE);limit=64*1024**2 if hard==resource.RLIM_INFINITY else min(hard,64*1024**2);resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
symbols=subprocess.check_output(['/usr/bin/nm','--defined-only',str(executable)],text=True).splitlines()
def sym(part):
    names=[x.split()[-1] for x in symbols if part in x and '.cold' not in x];assert len(names)==1;return names[0]
env=os.environ.copy();env.update(SIRIUS_GEOMETRY_CAPTURE=str(folder),SIRIUS_DISPATCH_SYMBOL=sym('VulkanDevice8DispatchE'),SIRIUS_LOAD_SYMBOL=sym('VulkanDevice10LoadKernelE'),SIRIUS_MODULE_INPUT=str(module),SIRIUS_FEATURE_LAYOUT=str(W/'feature-layout.json'),SIRIUS_FLOAT_LAYOUT=str(R/'bin/linux-gcc/lavapipe-native-fp32/float-layout.json'),SIRIUS_FEATURE_LOADER=str(Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1').resolve(strict=True)),SIRIUS_EXPECT_DENORM=str(int(arm=='candidate')),SIRIUS_ACTUAL_PROVIDER=str(driver))
script=folder/'capture.gdb';script.write_text('set pagination off\nset confirm off\nset breakpoint pending on\nset debuginfod enabled off\nset print thread-events off\nset disable-randomization off\nset environment LD_PRELOAD '+str(B/'optimizer_observer.so')+'\nset environment SIRIUS_OPT_IR_FOLDER '+str(folder)+'\nsource '+str(W/'capture_jit.py')+'\nsource '+str(W/'capture_features.py')+'\nsource '+str(W/'properties_capture.py')+'\nrun\n')
subprocess.run(['/usr/bin/gdb','-q','-batch','-nx','-x',str(script),'--args',str(executable),str(module),str(oracle/'input.bin'),str(oracle/'expected.bin'),str(folder/'actual.bin'),'--dispatch'],env=env,check=True)
assert not (folder/'capture-error.txt').exists()
code=json.loads((folder/'inferior-exit.json').read_text())['exit_code'];assert code in [0,1]
bindings=json.loads((folder/'bindings.json').read_text());assert len(bindings['dispatches'])==1 and bindings['dispatches'][0]['completed']
calls=[x for x in bindings['entry_calls'] if x['dispatch_id']==0];assert calls and len({x['worker_thread'] for x in calls})==len(calls)
assert all(x['actual_entry_pc']==x['function']['address'] for x in calls)
loaded=json.loads((folder/'actual-loaded-module.json').read_text());assert loaded['actual_span_matches_frozen_input'] and loaded['sha256']==hashlib.sha256(module.read_bytes()).hexdigest()
events=[x for x in bindings['jit_events'] if x.get('file') in {c['file'] for c in calls}];assert events and all(str(driver) in x['maps_at_first_call'] and 'libLLVM.so.21.1' in x['maps_at_first_call'] for x in events)
spec=importlib.util.spec_from_file_location('actual_features',W/'read_features.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader);features=reader.inspect(folder)
actual=(folder/'actual.bin').read_bytes();expected=(oracle/'expected.bin').read_bytes();assert len(actual)==len(expected)==522804*4
a=struct.unpack('<522804I',actual);e=struct.unpack('<522804I',expected);mismatches=[i for i in range(len(a)) if a[i]!=e[i]];untouched=sum(a[i]==(~e[i]&0xffffffff) for i in range(len(a)))
assert (code==0)==(not mismatches and untouched==0)
boundaries=json.loads((folder/'compute-boundaries.json').read_text());assert boundaries and all(x['returned'] for x in boundaries)
restored=all(x['mxcsr_before']==x['mxcsr_after'] for x in boundaries)
native_states=[{'worker':x['worker_thread'],'mxcsr':x['actual_mxcsr']} for x in calls]
if arm=='candidate':
    assert not mismatches and untouched==0 and restored
    assert all((x['mxcsr']&0xe040)==0 and (x['mxcsr']&0x1f80)==0x1f80 for x in native_states)
else:
    # Same exact baseline used by the later performance comparison must
    # reproduce the independently justified normal-domain rounding failure.
    assert poison and a[3*1508]==0x3f000001 and e[3*1508]==0x3f000000
    assert all(((x['mxcsr']>>13)&3)==2 for x in native_states)
result={'diagnostic_completed':True,'numerical_pass':not mismatches and untouched==0,'arm':arm,'caller_poisoned':poison,'mismatch_words':len(mismatches),'untouched_complement_words':untouched,'original_records':174268,'original_output_words':522804,'actual_module':loaded,'actual_features':features['features'],'worker_states':native_states,'callback_boundaries':boundaries,'all_observed_callback_states_restored':restored,'module_sha256':hashlib.sha256(module.read_bytes()).hexdigest(),'actual_provider_sha256':hashlib.sha256(driver.read_bytes()).hexdigest(),'scope':'Actual supported controls and actual implementation. Candidate additionally requests new genuinely implemented DenormPreserve32. Host modification is confined to the declared poisoned probe; observer is passive. No retained integration, benefit, scientific or release acceptance.'}
(folder/'numerical-and-boundary-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['arm','caller_poisoned','diagnostic_completed','numerical_pass','mismatch_words','all_observed_callback_states_restored']}))
