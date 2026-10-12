"""Compile the actual arithmetic witness with the production control recipe."""
from pathlib import Path
import importlib.util
import resource
import subprocess
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
work=Path(__file__).resolve().parent
build=root/'bin/linux-gcc/portable-execution-attribution/limb-screen'
assert len(sys.argv)==3
producer,mode=sys.argv[1:]
assert producer in ('baseline','candidate')
assert mode in ('integer','normal','normal_transport','normal_endpoint')
include=build/'baseline-kernels' if producer=='baseline' else root/'src/sirius/kernels'
stem=build/('arithmetic-'+producer+'-'+mode)
tools=root/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
raw=stem.with_suffix('.raw.spv')
optimized=stem.with_suffix('.optimized.spv')
assembly=stem.with_suffix('.spvasm')
output=stem.with_suffix('.spv')
assert not any(p.exists() for p in (raw,optimized,assembly,output))
definitions=['-DSIRIUS_RETAINED_PORTABLE=1']
if mode!='integer': definitions.append('-DSIRIUS_RETAINED_NORMAL_SUM32=1')
if mode=='normal_transport': definitions.append('-DSIRIUS_RETAINED_PARALLEL_TRANSPORT=1')
if mode=='normal_endpoint': definitions.append('-DSIRIUS_RETAINED_PARALLEL_ENDPOINT=1')
subprocess.run(['/opt/slang/bin/slangc',str(work/'limb_arithmetic_probe.slang'),*definitions,
    '-I',str(include),'-O0','-target','spirv','-profile','spirv_1_5',
    '-entry','ComputeMain','-stage','compute','-o',str(raw)],check=True)
subprocess.run(['/usr/bin/spirv-opt','--target-env=vulkan1.2','--ssa-rewrite',
    '--eliminate-dead-code-aggressive','--preserve-bindings','--preserve-interface',
    str(raw),'-o',str(optimized)],check=True)
subprocess.run([str(tools/'spirv-dis'),str(optimized),'-o',str(assembly)],check=True)
text=assembly.read_text()
if mode=='integer':
    assert 'OpTypeFloat' not in text
    output.write_bytes(optimized.read_bytes())
else:
    spec=importlib.util.spec_from_file_location('retained_controls',root/'scripts/build-retained-kernels.py')
    controls=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(controls)
    text=controls.portable_normal_sum_controls(text)
    assembly.write_text(text)
    subprocess.run([str(tools/'spirv-as'),'--target-env','spv1.5',str(assembly),'-o',str(output)],check=True)
subprocess.run([str(tools/'spirv-val'),'--target-env','vulkan1.2',str(output)],check=True)
print({'producer':producer,'mode':mode,'complete_module_valid':True,
       'normal_controls':'original production checker' if mode!='integer' else 'integer only'})
