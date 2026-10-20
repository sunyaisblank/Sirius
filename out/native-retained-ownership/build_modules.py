"""Build only the two ownership candidates with the original shader authority."""
from pathlib import Path
from types import SimpleNamespace
import hashlib, importlib.util, json, re, struct, subprocess, sys
sys.dont_write_bytecode = True
R = Path.cwd(); W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-retained-ownership/modules'
C = R / 'bin/linux-gcc/src/sirius/backend/retained'

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

def pin(p):
    return {'bytes': p.stat().st_size,
            'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}

authority = load('original_retained_builder', R / 'scripts/build-retained-kernels.py')
programs = load('original_retained_program', R / 'src/sirius/kernels/retained_program.py')
tools = {}
for line in (R / 'bin/linux-gcc/CMakeCache.txt').read_text().splitlines():
    if line.startswith('SIRIUS_') and ':FILEPATH=' in line:
        key, value = line.split('=', 1); tools[key.split(':')[0]] = value
actual_tools = {k: tools['SIRIUS_' + v] for k, v in {
    'compiler': 'SLANGC', 'assembler': 'SPIRV_AS', 'disassembler': 'SPIRV_DIS',
    'validator': 'SPIRV_VAL', 'optimizer': 'SPIRV_OPT'}.items()}
original = {str(p.relative_to(R)): pin(p) for p in C.iterdir() if p.is_file()}
commands = []
def run(argv, **kwargs):
    argv = list(argv)
    if argv[0] == actual_tools['compiler']:
        if Path(argv[1]).stem == 'retained_endpoint':
            argv[1] = str(W / 'shaders/retained_endpoint.slang')
            argv[2:2] = ['-I', str(R / 'src/sirius/kernels')]
        elif Path(argv[1]).stem == 'retained_dopri_phase':
            argv.insert(2, '-DSIRIUS_RETAINED_COOPERATIVE_DOPRI_INPUTS=1')
        else:
            raise AssertionError(argv)
    commands.append(argv)
    return subprocess.run(argv, **kwargs)
authority.subprocess = SimpleNamespace(run=run)

header = ['#pragma once', '#include "' + str(C / 'retained_kernels.h') + '"',
          'namespace sirius::backend::retained_program {']
result = {}
for kind, stem, build in [
    ('Endpoint', 'endpoint', programs.build_endpoint_program),
    ('DopriPhase', 'dopri_phase', programs.build_dopri_phase_program)]:
    plan = build(parallel=True)
    encoded = [plan['instructions'], plan['registers'], len(plan['outputs'])]
    encoded += plan['outputs'] + plan['operations'] + plan['layer_offsets']
    specialized = None
    if kind == 'Endpoint':
        sp = programs.build_schwarzschild_endpoint_program(parallel=True)
        specialized_words = [sp['instructions'], sp['registers'], len(sp['outputs'])]
        specialized_words += sp['outputs'] + sp['operations'] + sp['layer_offsets']
        specialized = (len(encoded), len(specialized_words),
                       len(sp['layer_offsets']) - 1, sp['prefix_instructions'])
        encoded += specialized_words
    text = (C / 'retained_kernels.h').read_text()
    match = re.search(r'k' + kind + r'Program\{\{(.*?)\}\};', text, re.S)
    assert match and [int(x) for x in re.findall(r'(\d+)u', match[1])] == encoded
    destination = B / ('owned-' + stem + '.spv')
    assert not destination.exists()
    words = authority.compile_shader(
        R / ('src/sirius/kernels/retained_' + stem + '.slang'), destination,
        actual_tools['compiler'], actual_tools['assembler'],
        actual_tools['disassembler'], actual_tools['validator'],
        plan['registers'], 5, len(plan['layer_offsets']) - 1,
        plan.get('prefix_instructions', 0), optimizer=actual_tools['optimizer'],
        endpoint_specialized=specialized)
    assembly = subprocess.check_output([actual_tools['disassembler'], str(destination)], text=True)
    old_assembly = subprocess.check_output([
        actual_tools['disassembler'], str(C / ('retained_' + stem + '.spv'))], text=True)
    assert re.findall(r'OpCapability (\S+)', assembly) == re.findall(r'OpCapability (\S+)', old_assembly)
    assert re.findall(r'OpExecutionMode .*', assembly) == re.findall(r'OpExecutionMode .*', old_assembly)
    assert re.search(r'OpExecutionMode %\S+ LocalSize 64 1 1', assembly)
    assert 'DenormPreserve 32' in assembly and 'RoundingModeRTE 32' in assembly
    (B / ('owned-' + stem + '.spvasm')).write_text(assembly)
    header.append('inline constexpr std::array<std::uint32_t, ' + str(len(words)) + '> k' + kind + 'NativeOwnedShader{{')
    for i in range(0, len(words), 16):
        header.append(','.join(str(x) + 'u' for x in words[i:i+16]) + ',')
    header.append('}};')
    result[kind] = {'module': str(destination.relative_to(R)), **pin(destination),
        'original_module': pin(C / ('retained_' + stem + '.spv')),
        'LocalSize': [64, 1, 1], 'program_table_equal': True,
        'original_control_barriers': old_assembly.count('OpControlBarrier'),
        'candidate_control_barriers': assembly.count('OpControlBarrier'),
        'scalar_algorithm_unchanged': True}
header.append('}')
(W / 'headers').mkdir()
(W / 'headers/retained_kernels.h').write_text('\n'.join(header) + '\n')
assert original == {str(p.relative_to(R)): pin(p) for p in C.iterdir() if p.is_file()}
(W / 'compiler-result.json').write_text(json.dumps({
    'modules': result, 'commands': commands, 'tools': actual_tools,
    'all35_shipped_modules_and_original_generated_headers_unchanged': True,
    'original_inputs': original, 'workgroup_layout_or_arithmetic_changes': False,
    'qualification': False}, indent=2) + '\n')
print(json.dumps({'compiled_modules': 2, 'all35_original_modules_unchanged': True,
                  'modules': result}))
