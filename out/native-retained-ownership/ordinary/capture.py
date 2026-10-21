from pathlib import Path
import hashlib,importlib.util,json,os,resource,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;kind=sys.argv[1];folder=Path(sys.argv[2]);B=R/'bin/linux-gcc/native-retained-ownership';driver=B/'providers/candidate/libvulkan_lvp.so';exe=B/'short_input_jit_probe';module=R/json.loads((W.parent/'compiler-result.json').read_text())['modules'][kind]['module']
assert json.loads((W.parent/'reused/prior-proof/lavapipe-ieee32-validation/candidate-default/numerical-and-boundary-result.json').read_text())['numerical_pass']
assert json.loads((W.parent/'reused/prior-proof/lavapipe-ieee32-validation/candidate-poison/numerical-and-boundary-result.json').read_text())['numerical_pass']
resource.setrlimit(resource.RLIMIT_CORE,(0,0));_,hard=resource.getrlimit(resource.RLIMIT_FSIZE);limit=64*1024**2 if hard==resource.RLIM_INFINITY else min(hard,64*1024**2);resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
symbols=subprocess.check_output(['/usr/bin/nm','--defined-only',str(exe)],text=True).splitlines()
def sym(part):
    found=[x.split()[-1] for x in symbols if part in x and '.cold' not in x];assert len(found)==1;return found[0]
env=os.environ.copy();env.update(SIRIUS_GEOMETRY_CAPTURE=str(folder),SIRIUS_DISPATCH_SYMBOL=sym('VulkanDevice8DispatchE'),SIRIUS_LOAD_SYMBOL=sym('VulkanDevice10LoadKernelE'),SIRIUS_MODULE_INPUT=str(module),SIRIUS_FEATURE_LAYOUT=str(W/'feature-layout.json'),SIRIUS_FLOAT_LAYOUT=str(W/'float-layout.json'),SIRIUS_FEATURE_LOADER=str(Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1').resolve(strict=True)),SIRIUS_EXPECT_DENORM='1',SIRIUS_ACTUAL_PROVIDER=str(driver))
script=folder/'capture.gdb';script.write_text('set pagination off\nset confirm off\nset breakpoint pending on\nset debuginfod enabled off\nset print thread-events off\nset disable-randomization off\nset environment LD_PRELOAD '+str(B/'optimizer_observer.so')+'\nset environment SIRIUS_OPT_IR_FOLDER '+str(folder)+'\nsource '+str(W/'capture_jit.py')+'\nsource '+str(W/'capture_features.py')+'\nsource '+str(W/'properties_capture.py')+'\nrun\n')
subprocess.run(['/usr/bin/gdb','-q','-batch','-nx','-x',str(script),'--args',str(exe),str(module)],env=env,check=True)
assert not (folder/'capture-error.txt').exists() and json.loads((folder/'inferior-exit.json').read_text())['exit_code']==0
bindings=json.loads((folder/'bindings.json').read_text());assert len(bindings['dispatches'])==1 and bindings['dispatches'][0]['completed']
entry,=[x for x in bindings['entry_calls'] if x['dispatch_id']==0]
event,=[x for x in bindings['jit_events'] if x.get('file')==entry['file']]
assert str(driver) in event['maps_at_first_call'] and 'libLLVM.so.21.1' in event['maps_at_first_call']
lookups=[x for x in bindings['pointer_lookups'] if entry['file'] in x['matched_files']];generations={x['module_generation'] for x in lookups};assert len(generations)==1
export,=[x for x in bindings['ir_exports'] if x['module_generation']==next(iter(generations))]
assert export['optimized_module_observed'] and all(bindings['optimizer_passes'][i]['completed'] and bindings['optimizer_passes'][i]['error_pointer']==0 for i in export['successful_optimizer_pass_ids'])
loaded=json.loads((folder/'actual-loaded-module.json').read_text());assert loaded['actual_span_matches_frozen_input'] and loaded['sha256']==hashlib.sha256(module.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('feature_reader',W/'read_features.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader);features=reader.inspect(folder)
boundary=json.loads((folder/'compute-boundaries.json').read_text());assert boundary and all(x['returned'] and x['mxcsr_before']==x['mxcsr_after'] for x in boundary)
(folder/'capture-verification.json').write_text(json.dumps({'passed':True,'actual_entry':entry,'ordinary_ir_export':export,'actual_complete_module':loaded,'actual_logical_features':features['features'],'actual_provider_sha256':hashlib.sha256(driver.read_bytes()).hexdigest(),'callback_restoration_observed':True,'scope':'Actual complete ownership candidate module and ordinary called entry,4B guard/canary dispatch. No retained numerical path or benefit/frame/release claim.'},indent=2)+'\n');print(json.dumps({'passed':True,'module_bytes':loaded['bytes'],'total_JIT_function_bytes':sum(x['bytes'] for x in event['functions'])}))
