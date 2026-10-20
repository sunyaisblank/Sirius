from pathlib import Path
import hashlib,io,json,shutil,subprocess,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-ieee32-staged'
assert json.loads((R/'out/lavapipe-ieee32/build/owner.json').read_text())['passed']
c=B/'candidate';c.mkdir();source=R/'bin/linux-gcc/lavapipe-ieee32/src/gallium/targets/lavapipe/libvulkan_lvp.so';shutil.copyfile(source,c/'libvulkan_lvp.so')
m=json.loads((B/'baseline/lvp_icd.json').read_text());m['ICD']['library_path']=str(c/'libvulkan_lvp.so');(c/'lvp_icd.json').write_text(json.dumps(m,indent=2)+'\n')
prefix='attestations/qualification/8c2be67d/lavapipe-vector-workgroup-rejected/'
with zipfile.ZipFile(io.BytesIO(subprocess.check_output(['git','show','43389d91:'+prefix+'providers-and-probes.zip']))) as archive:
    (B/'optimizer_observer.so').write_bytes(archive.read('bin/linux-gcc/lavapipe-vector-inspect/optimizer_observer.so'))
with zipfile.ZipFile(io.BytesIO(subprocess.check_output(['git','show','43389d91:'+prefix+'diagnostics-and-trial.zip']))) as archive:
    (W/'optimizer_observer.c').write_bytes(archive.read('out/lavapipe-vector-inspect/optimizer_observer.c'))
native=R/'bin/linux-gcc/lavapipe-native-fp32';s=(native/'native.spvasm').read_text()
assert 'DenormPreserve' not in s;s=s.replace('OpCapability Shader','OpCapability Shader\n               OpCapability DenormPreserve',1)
import re
entry=re.search(r'OpEntryPoint GLCompute (%\S+)',s)[1];line=re.search(r'^.*OpExecutionMode '+re.escape(entry)+r' LocalSize.*$',s,re.M)[0]
s=s.replace(line,line+'\n               OpExecutionMode '+entry+' DenormPreserve 32',1)
(B/'native-preserve.spvasm').write_text(s);tools=json.loads((R/'out/lavapipe-native-fp32/retained-endpoint-arguments.json').read_text())['tools']
subprocess.run([tools['assembler'],'--target-env','spv1.5',str(B/'native-preserve.spvasm'),'-o',str(B/'native-preserve.spv')],check=True);subprocess.run([tools['validator'],'--target-env','vulkan1.2',str(B/'native-preserve.spv')],check=True)
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
providers={arm:{'library':str(B/arm/'libvulkan_lvp.so'),'ICD':str(B/arm/'lvp_icd.json'),'library_identity':pin(B/arm/'libvulkan_lvp.so'),'ICD_identity':pin(B/arm/'lvp_icd.json')} for arm in ['baseline','candidate']}
(W/'providers.json').write_text(json.dumps({'providers':providers,'observer':pin(B/'optimizer_observer.so'),'candidate_equals_actual_built_library':pin(source)==pin(c/'libvulkan_lvp.so'),'baseline_provider_reused_from_pushed43389d91':True,'default_module':pin(native/'native.spv'),'preserve_module':pin(B/'native-preserve.spv'),'instruction_delta':'Add only DenormPreserve capability/execution-mode32 to the candidate primitive module; operations and original expected words unchanged.','system_provider_or_toolchain_installations':0},indent=2)+'\n')
for name in ['run_owned.py','owned_guard.py']:
    shutil.copyfile(R/'out/lavapipe-native-fp32'/name,W/name)
p=W/'run_owned.py';s=p.read_text();old="'VK_DRIVER_FILES':'/usr/share/vulkan/icd.d/lvp_icd.json'";assert s.count(old)==1;s=s.replace(old,"'VK_DRIVER_FILES':action['runtime_icd']");p.write_text(s)
