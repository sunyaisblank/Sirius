from pathlib import Path
import hashlib,json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/lavapipe-native-fp32'
t=json.loads((W/'retained-endpoint-arguments.json').read_text())['tools']
subprocess.run(['/usr/bin/g++-14','-std=c++2c','-Wall','-Wextra','-Werror',str(W/'layout.cpp'),'-o',str(B/'layout')],check=True)
layout=json.loads(subprocess.check_output([str(B/'layout')],text=True));assert layout['first_field_offset']==16 and layout['last_field_offset']==80 and layout['bytes']==88
(B/'float-layout.json').write_text(json.dumps(layout,indent=2)+'\n')
subprocess.run([t['compiler'],str(W/'corpus.slang'),'-I',str(W),'-O0','-target','spirv','-profile','spirv_1_5','-entry','main','-stage','compute','-o',str(B/'native.compiler.spv')],check=True)
subprocess.run([t['optimizer'],'--target-env=vulkan1.2','--ssa-rewrite','--eliminate-dead-code-aggressive','--preserve-bindings','--preserve-interface',str(B/'native.compiler.spv'),'-o',str(B/'native.optimized.spv')],check=True)
subprocess.run([t['disassembler'],str(B/'native.optimized.spv'),'-o',str(B/'native.spvasm')],check=True)
p=B/'native.spvasm';s=p.read_text();assert re.findall(r'OpCapability (\S+)',s)==['Shader'] and set(re.findall(r'OpTypeFloat (\d+)',s))=={'32'} and set(re.findall(r'OpTypeInt (\d+) [01]',s))=={'32'}
entry=re.search(r'OpEntryPoint GLCompute (%\S+)',s)[1]
s=s.replace('OpCapability Shader','OpCapability Shader\n               OpCapability RoundingModeRTE\n               OpCapability SignedZeroInfNanPreserve',1)
local=re.search(r'^.*OpExecutionMode '+re.escape(entry)+r' LocalSize.*$',s,re.M)[0]
s=s.replace(local,local+'\n               OpExecutionMode '+entry+' RoundingModeRTE 32\n               OpExecutionMode '+entry+' SignedZeroInfNanPreserve 32',1)
assert 'DenormPreserve' not in s and 'DenormFlushToZero' not in s
p.write_text(s)
subprocess.run([t['assembler'],'--target-env','spv1.5',str(p),'-o',str(B/'native.spv')],check=True)
subprocess.run([t['validator'],'--target-env','vulkan1.2',str(B/'native.spv')],check=True)
(W/'compiler-result.json').write_text(json.dumps({'passed':True,'module_sha256':hashlib.sha256((B/'native.spv').read_bytes()).hexdigest(),'only_native_binary32':True,'integer_product_residual_unchanged':True,'denormal_execution_mode_requested':False,'required_modes':['RoundingModeRTE32','SignedZeroInfNanPreserve32'],'device_calls':0,'not_a_numerical_or_feature_qualification':True},indent=2)+'\n')
