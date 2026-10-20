"""One observed default-native scalar dispatch; mismatch is evidence, not a pass."""
from pathlib import Path
import collections,hashlib,importlib.util,json,os,re,resource,struct,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-native-fp32';folder=Path(sys.argv[1]).resolve(strict=True)
assert folder.parent==W
assert json.loads((W/'build/owner.json').read_text())['passed']
module=B/'native.spv';executable=B/'gpu_probe'
observer=R/'attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/optimizer_observer.so'
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
_,hard=resource.getrlimit(resource.RLIMIT_FSIZE);limit=64*1024**2 if hard==resource.RLIM_INFINITY else min(hard,64*1024**2);resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
symbols=subprocess.check_output(['/usr/bin/nm','--defined-only',str(executable)],text=True).splitlines()
def sym(part):
    matches=[x.split()[-1] for x in symbols if part in x and '.cold' not in x];assert len(matches)==1,matches;return matches[0]
env=os.environ.copy();env.update(SIRIUS_GEOMETRY_CAPTURE=str(folder),SIRIUS_DISPATCH_SYMBOL=sym('VulkanDevice8DispatchE'),SIRIUS_LOAD_SYMBOL=sym('VulkanDevice10LoadKernelE'),SIRIUS_MODULE_INPUT=str(module),SIRIUS_FEATURE_LAYOUT=str(W/'feature-layout.json'),SIRIUS_FLOAT_LAYOUT=str(B/'float-layout.json'),SIRIUS_FEATURE_LOADER=str(Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1').resolve(strict=True)))
script=folder/'capture.gdb';script.write_text('set pagination off\nset confirm off\nset breakpoint pending on\nset debuginfod enabled off\nset print thread-events off\nset disable-randomization off\nset environment LD_PRELOAD '+str(observer)+'\nset environment SIRIUS_OPT_IR_FOLDER '+str(folder)+'\nsource '+str(W/'capture_jit.py')+'\nsource '+str(W/'capture_features.py')+'\nsource '+str(W/'properties_capture.py')+'\nrun\n')
subprocess.run(['/usr/bin/gdb','-q','-batch','-nx','-x',str(script),'--args',str(executable),str(module),str(W/'input.bin'),str(W/'expected.bin'),str(folder/'actual.bin'),'--dispatch'],env=env,check=True)
assert not (folder/'capture-error.txt').exists()
exit_code=json.loads((folder/'inferior-exit.json').read_text())['exit_code'];assert exit_code in [0,1]
bindings=json.loads((folder/'bindings.json').read_text());assert len(bindings['dispatches'])==1 and bindings['dispatches'][0]['completed']
calls=[x for x in bindings['entry_calls'] if x['dispatch_id']==0];assert calls and len({x['worker_thread'] for x in calls})==len(calls)
assert all(x['function']['address']==x['actual_entry_pc'] for x in calls)
loaded=json.loads((folder/'actual-loaded-module.json').read_text());assert loaded['actual_span_matches_frozen_input'] and loaded['sha256']==hashlib.sha256(module.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('actual_feature_reader',W/'read_features.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader);features=reader.inspect(folder)
properties=json.loads((folder/'float-controls-properties.json').read_text());assert properties and all(x['values']['shaderDenormPreserveFloat32']==0 and x['values']['shaderRoundingModeRTEFloat32']==x['values']['shaderSignedZeroInfNanPreserveFloat32']==1 for x in properties)
expected=(W/'expected.bin').read_bytes();actual=(folder/'actual.bin').read_bytes();assert len(actual)==len(expected)==522804*4
e=struct.unpack('<522804I',expected);a=struct.unpack('<522804I',actual);inputs=struct.unpack('<'+str((W/'input.bin').stat().st_size//4)+'I',(W/'input.bin').read_bytes())
different=[i for i in range(len(e)) if a[i]!=e[i]];untouched=sum(a[i]==(~e[i]&0xffffffff) for i in range(len(e)))
assert (exit_code==0)==(not different and untouched==0)
counts=collections.Counter(inputs[1+3*(i//3)] for i in different)
examples=[{'row':i//3,'column':i%3,'operation':inputs[1+3*(i//3)],'a':hex(inputs[2+3*(i//3)]),'b':hex(inputs[3+3*(i//3)]),'expected':hex(e[i]),'actual':hex(a[i])} for i in different[:20]]
states=[{'worker':x['worker_thread'],'mxcsr':hex(x['actual_mxcsr']),'FTZ':bool(x['actual_mxcsr']&0x8000),'DAZ':bool(x['actual_mxcsr']&0x40),'rounding_control':(x['actual_mxcsr']>>13)&3} for x in calls]
result={'diagnostic_completed':True,'original_records':174268,'complete_output_words':522804,'native_arithmetic_matches_original':not different and untouched==0,'mismatch_words':len(different),'untouched_complement_words':untouched,'mismatches_by_operation':dict(counts),'examples':examples,'actual_workers':states,'public_feature_values':features['features'],'actual_control_properties':properties,'complete_module':loaded,'scope':'One real stock-provider default-denormal dispatch. Actual supported RTE32/SignedZeroInfNanPreserve32 only; no DenormPreserve mode or property/register override. This observation is not a feature guarantee, retained integration, benefit, frame or release qualification.'}
(folder/'native-default-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['diagnostic_completed','native_arithmetic_matches_original','mismatch_words','mismatches_by_operation','actual_workers','scope']}))
