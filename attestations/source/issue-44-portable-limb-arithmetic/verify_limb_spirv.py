"""Verify all accumulation specializations in the already compiled products."""
from pathlib import Path
import hashlib
import json
import re
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/portable-execution-attribution/limb-screen'
tools = root / 'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
paths = {'baseline':root / 'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv',
         'candidate':stage / 'retained_endpoint_portable_normal_sum.spv'}
reports = []
for producer,path in paths.items():
    subprocess.run([str(tools/'spirv-val'),'--target-env','vulkan1.2',str(path)],check=True)
    assembly = stage / (producer+'-verified.spvasm')
    assert not assembly.exists()
    subprocess.run([str(tools/'spirv-dis'),str(path),'-o',str(assembly)],check=True)
    text=assembly.read_text()
    bodies=re.findall(r'(%\S+) = OpFunction .*?OpFunctionEnd',text,re.S)
    accumulations=[(match.group(1),match.group(0)) for match in
                   re.finditer(r'(%\S+) = OpFunction .*?OpFunctionEnd',text,re.S)
                   if re.fullmatch(r'%RTAccumulateTerm(?:_\d+)?',match.group(1))]
    assert len(accumulations)==2, (producer,'sealed module specializations')
    details=[]
    for name,body in accumulations:
        calls=re.findall(r'= OpFunctionCall %\S+ (%\S+)',body)
        sums=[value for value in calls if value.startswith('%RetainedSum32')]
        radii=[value for value in calls if value.startswith('%RPUpAdd')]
        bounds=[value for value in calls if value.startswith('%RPBounded')]
        assert len(set(sums))==1
        assert len(sums)==(1 if producer=='baseline' else 3)
        assert body.count('OpLoopMerge')==(1 if producer=='baseline' else 0)
        assert len(radii)==1 and len(bounds)==5
        details.append({'function':name,'loop_merges':body.count('OpLoopMerge'),
                        'sum_calls':sums,'radius_calls':radii,'bound_calls':bounds,
                        'complete_calls':calls})
    entry=re.search(r'OpEntryPoint GLCompute (%\S+)',text).group(1)
    assert f'OpExecutionMode {entry} LocalSize 16 1 1' in text
    assert f'OpExecutionMode {entry} RoundingModeRTE 32' in text
    assert set(re.findall(r'OpTypeInt (\d+) [01]',text))=={'32'}
    assert re.findall(r'OpTypeFloat (\d+)',text)==['32']
    assert re.findall(r'OpCapability (\S+)',text)==['Shader','RoundingModeRTE']
    reports.append({'producer':producer,'module_bytes':path.stat().st_size,
                    'module_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                    'specializations':details,'whole_module_loop_merges':text.count('OpLoopMerge'),
                    'whole_module_functions':len(bodies),
                    'control_barriers':len(re.findall(r'\bOpControlBarrier\b',text)),
                    'memory_barriers':len(re.findall(r'\bOpMemoryBarrier\b',text))})
assert reports[0]['control_barriers']==reports[1]['control_barriers']
assert reports[0]['memory_barriers']==reports[1]['memory_barriers']
assert reports[1]['whole_module_loop_merges']==reports[0]['whole_module_loop_merges']-2
assert reports[1]['whole_module_functions']==reports[0]['whole_module_functions']
result={'passed':True,'compiler_rerun':False,'reports':reports,
        'scope':'Exact retained normal Endpoint products validate. Both accumulation specializations have three explicit sum calls, no limb loop, unchanged single radius/five bound calls and original modes/layout/barriers. No machine/private-frame/numerical-runtime/performance/full-products/qualification claim.'}
destination=work/'limb-screen-spirv-verification.json'
assert not destination.exists()
destination.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
