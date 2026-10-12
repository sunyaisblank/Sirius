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
stage = root / 'bin/linux-gcc/portable-indexed-limb-accumulation'
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
assert output.is_file()
result = {'passed': True, 'module_path': str(output.relative_to(root)),
          'module_bytes': output.stat().st_size,
          'module_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
          'program_registers': general['registers'], 'program_layers': len(general['layer_offsets'])-1,
          'program_prefix': general.get('prefix_instructions',0), 'specialized_layout': layout,
          'scope': 'One staged original normal Endpoint with unchanged production compiler/SSA/controls/validation recipe. No dispatch, numerical, machine/private-frame, performance or qualification claim.'}
(work/'compiler-product.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
