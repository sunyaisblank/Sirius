"""Compile one staged normal Endpoint with the original production recipe."""
from pathlib import Path
import hashlib
import importlib.util
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
source = stage / 'kernels'
tools = root / 'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value
controls = load('production_retained_controls', root / 'scripts/build-retained-kernels.py')
programs = load('production_retained_program', source / 'retained_program.py')
general = programs.build_endpoint_program(parallel=True)
special = programs.build_schwarzschild_endpoint_program(parallel=True)
assert general['registers'] == special['registers'] and len(general['outputs']) == len(special['outputs']) == 100
assert special['prefix_instructions'] in special['layer_offsets']
encode = lambda p: [p['instructions'],p['registers'],100]+p['outputs']+p['operations']+p['layer_offsets']
layout = (len(encode(general)), len(encode(special)), len(special['layer_offsets'])-1, special['prefix_instructions'])
output = stage / 'retained_endpoint_portable_normal_sum.spv'
assert not output.exists()
code = controls.compile_shader(source / 'retained_endpoint.slang', output,
    '/opt/slang/bin/slangc', str(tools / 'spirv-as'), str(tools / 'spirv-dis'), str(tools / 'spirv-val'),
    general['registers'], 5, len(general['layer_offsets'])-1, general.get('prefix_instructions',0),
    portable=True, optimizer='/usr/bin/spirv-opt', normal_sum32=True, endpoint_specialized=layout)
baseline = root / 'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv'
reports = []
for producer, path in [('baseline', baseline), ('candidate', output)]:
    assembly = stage / (producer + '.spvasm')
    assert not assembly.exists()
    subprocess.run([str(tools / 'spirv-dis'), str(path), '-o', str(assembly)], check=True)
    text = assembly.read_text()
    bodies = re.findall(r'%\S+ = OpFunction .*?OpFunctionEnd', text, re.S)
    accumulations = [body for body in bodies if re.match(r'%RTAccumulateTerm(?:_\d+)? = OpFunction ',body)]
    assert len(accumulations)==1
    body = accumulations[0]
    calls = re.findall(r'= OpFunctionCall %\S+ (%\S+)', body)
    reports.append({'producer':producer,'module_bytes':path.stat().st_size,
                    'module_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                    'RTAccumulateTerm_loop_merges':body.count('OpLoopMerge'),
                    'RTAccumulateTerm_calls':calls,
                    'whole_module_loop_merges':text.count('OpLoopMerge'),
                    'whole_module_functions':len(bodies)})
assert reports[0]['RTAccumulateTerm_loop_merges']==1 and reports[1]['RTAccumulateTerm_loop_merges']==0
assert sum(name.startswith('%RetainedSum32') for name in reports[0]['RTAccumulateTerm_calls'])==1
assert sum(name.startswith('%RetainedSum32') for name in reports[1]['RTAccumulateTerm_calls'])==3
result={'passed':True,'scope':'One actual normal Endpoint compiler feasibility screen only; SPIR-V validates with original controls/layout/recipe. No machine-effect, private-frame, arithmetic/runtime/performance or full-product/release qualification claim.',
        'program_registers':general['registers'],'program_layers':len(general['layer_offsets'])-1,
        'program_prefix':general.get('prefix_instructions',0),'specialized_layout':layout,
        'reports':reports}
destination=work/'limb-screen-compiler-verification.json'
assert not destination.exists()
destination.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
