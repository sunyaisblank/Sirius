"""Run the reused passive JIT/feature capture against one exact private provider."""
from pathlib import Path
import hashlib,importlib.util,json,os,resource,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;P=R/'out/lavapipe-vector-inspect'
folder=Path(sys.argv[1]).resolve(strict=True);arm=sys.argv[2]
assert arm in ['baseline','candidate'] and folder.parent==W
module=R/'bin/linux-gcc/src/sirius/backend/retained/retained_transport_portable_normal_sum.spv'
executable=R/'bin/linux-gcc/lavapipe-vector-inspect/short_input_jit_probe'
observer=R/'bin/linux-gcc/lavapipe-vector-inspect/optimizer_observer.so'
driver=R/f'bin/linux-gcc/lavapipe-vector-staged/{arm}/libvulkan_lvp.so'
icd=P/f'{arm}-icd.json'
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
_,hard=resource.getrlimit(resource.RLIMIT_FSIZE);limit=64*1024**2 if hard==resource.RLIM_INFINITY else min(hard,64*1024**2)
resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
symbols=subprocess.check_output(['/usr/bin/nm','--defined-only',str(executable)],text=True).splitlines()
def sym(part):
    matches=[line.split()[-1] for line in symbols if part in line and '.cold' not in line];assert len(matches)==1,matches
    return matches[0]
env=os.environ.copy();env.update(VK_DRIVER_FILES=str(icd),SIRIUS_VULKAN_DEVICE='0',MESA_SHADER_CACHE_DISABLE='true',
 SIRIUS_GEOMETRY_CAPTURE=str(folder),SIRIUS_DISPATCH_SYMBOL=sym('VulkanDevice8DispatchE'),
 SIRIUS_LOAD_SYMBOL=sym('VulkanDevice10LoadKernelE'),SIRIUS_MODULE_INPUT=str(module),
 SIRIUS_FEATURE_LAYOUT=str(P/'feature-layout.json'),SIRIUS_FEATURE_LOADER=str(Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1').resolve(strict=True)))
script=folder/'capture.gdb'
script.write_text('set pagination off\nset confirm off\nset breakpoint pending on\nset debuginfod enabled off\n'
 'set print thread-events off\nset disable-randomization off\nset environment LD_PRELOAD '+str(observer)+'\n'
 'set environment SIRIUS_OPT_IR_FOLDER '+str(folder)+'\nsource '+str(P/'capture_jit.py')+'\nsource '+str(P/'capture_features.py')+'\nrun\n')
(folder/'capture-inputs.json').write_text(json.dumps({'arm':arm,'module':str(module),'module_sha256':sha(module),
 'driver':str(driver),'driver_sha256':sha(driver),'ICD':str(icd),'ICD_sha256':sha(icd),'host_sha256':sha(executable),
 'observer_sha256':sha(observer),'runtime_environment':{k:env[k] for k in ['VK_DRIVER_FILES','SIRIUS_VULKAN_DEVICE','MESA_SHADER_CACHE_DISABLE']},
 'scope':'Unchanged complete Transport module, one-word input guard, one public factory dispatch; no retained arithmetic or timing verdict.'},indent=2)+'\n')
subprocess.run(['/usr/bin/gdb','-q','-batch','-nx','-x',str(script),'--args',str(executable),str(module)],env=env,check=True)
assert not (folder/'capture-error.txt').exists()
assert json.loads((folder/'inferior-exit.json').read_text())['exit_code']==0
bindings=json.loads((folder/'bindings.json').read_text())
assert len(bindings['dispatches'])==1 and bindings['dispatches'][0]['completed']
calls=[x for x in bindings['entry_calls'] if x['dispatch_id']==0];assert len(calls)==1
lookups=[x for x in bindings['pointer_lookups'] if calls[0]['file'] in x['matched_files']]
generations={x['module_generation'] for x in lookups};assert len(generations)==1
generation=generations.pop();exports=[x for x in bindings['ir_exports'] if x['module_generation']==generation]
assert len(exports)==1 and exports[0]['optimized_module_observed'];export=exports[0]
passes=[bindings['optimizer_passes'][i] for i in export['successful_optimizer_pass_ids']]
assert all(x['completed'] and x['error_pointer']==0 for x in passes)
assert [x['pipeline'] for x in passes]==['default<O0>','sroa,early-cse,simplifycfg,reassociate,mem2reg,instsimplify,instcombine<no-verify-fixpoint>']
loaded=json.loads((folder/'actual-loaded-module.json').read_text());assert loaded['actual_span_matches_frozen_input'] and loaded['sha256']==sha(module)
event,=[x for x in bindings['jit_events'] if x.get('file')==calls[0]['file']]
maps=event['maps_at_first_call'];assert str(driver) in maps
assert '/libLLVM.so.21.1' in maps and '/libvulkan.so.1.3.275' in maps
(folder/'capture-verification.json').write_text(json.dumps({'passed':True,'actual_loaded_module':loaded,'ordinary_ir_export':export,
 'actual_entry':calls[0],'module_generation':generation,'actual_provider_mapping':True,
 'scope':'Actual unchanged loaded module, completed short-input dispatch and called entry under matched private provider. No retained arithmetic or benefit qualification.'},indent=2)+'\n')
spec=importlib.util.spec_from_file_location('ordinary_feature_reader',P/'read_features.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
(folder/'feature-verification.json').write_text(json.dumps(reader.inspect(folder),indent=2)+'\n')
print(json.dumps({'passed':True,'arm':arm,'called_entry_bytes':calls[0]['function']['bytes'],'IR':export['path']}))
