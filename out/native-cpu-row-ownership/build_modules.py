"""Compile three private complete-row modules; preserve the original tables/ABI."""
from pathlib import Path
from types import SimpleNamespace
import hashlib, importlib.util, json, re, subprocess, sys
sys.dont_write_bytecode = True
R = Path.cwd(); W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-row-ownership/modules'
C = R / 'bin/linux-gcc/src/sirius/backend/retained'
def load(name, p):
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
def pin(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def write(p, d): p.write_text(json.dumps(d, indent=2) + '\n')
authority = load('original_builder', R / 'scripts/build-retained-kernels.py')
programs = load('original_programs', R / 'src/sirius/kernels/retained_program.py')
tools = {}
for line in (R / 'bin/linux-gcc/CMakeCache.txt').read_text().splitlines():
    if line.startswith('SIRIUS_') and ':FILEPATH=' in line:
        k, v = line.split('=', 1); tools[k.split(':')[0]] = v
tool = {k: tools['SIRIUS_' + v] for k, v in {'compiler': 'SLANGC', 'assembler': 'SPIRV_AS',
    'disassembler': 'SPIRV_DIS', 'validator': 'SPIRV_VAL', 'optimizer': 'SPIRV_OPT'}.items()}
original = {str(p.relative_to(R)): pin(p) for p in C.iterdir() if p.is_file()}
commands = []
def run(argv, **kw):
    argv = list(argv)
    if argv[0] == tool['compiler']:
        assert Path(argv[1]).stem in ['retained_transport', 'retained_endpoint', 'retained_dopri_phase']
        argv[1] = str(W / 'shaders' / Path(argv[1]).name)
        argv[2:2] = ['-DSIRIUS_RETAINED_CPU_ROW_OWNER=1', '-I', str(W / 'shaders'),
                    '-I', str(R / 'src/sirius/kernels')]
    commands.append(argv); return subprocess.run(argv, **kw)
authority.subprocess = SimpleNamespace(run=run)
header = ['#pragma once', '#include "' + str(C / 'retained_kernels.h') + '"',
          'namespace sirius::backend::retained_program {']
rows = {}
for kind, stem, build in [('Transport', 'transport', programs.build_transport_program),
                          ('Endpoint', 'endpoint', programs.build_endpoint_program),
                          ('DopriPhase', 'dopri_phase', programs.build_dopri_phase_program)]:
    plan = build(parallel=True)
    encoded = [plan['instructions'], plan['registers'], len(plan['outputs'])]
    encoded += plan['outputs'] + plan['operations'] + plan['layer_offsets']
    kwargs = {}
    if kind in ['Transport', 'Endpoint']:
        sp = (programs.build_schwarzschild_transport_program if kind == 'Transport'
              else programs.build_schwarzschild_endpoint_program)(parallel=True)
        extra = [sp['instructions'], sp['registers'], len(sp['outputs'])]
        extra += sp['outputs'] + sp['operations'] + sp['layer_offsets']
        frame = (len(encoded), len(extra), len(sp['layer_offsets']) - 1)
        if kind == 'Endpoint': frame += (sp['prefix_instructions'],)
        kwargs['transport_specialized' if kind == 'Transport' else 'endpoint_specialized'] = frame
        encoded += extra
    match = re.search(r'k' + kind + r'Program\{\{(.*?)\}\};', (C / 'retained_kernels.h').read_text(), re.S)
    assert match and [int(v) for v in re.findall(r'(\d+)u', match[1])] == encoded
    destination = B / ('owned-' + stem + '.spv'); assert not destination.exists()
    words = authority.compile_shader(R / ('src/sirius/kernels/retained_' + stem + '.slang'),
        destination, tool['compiler'], tool['assembler'], tool['disassembler'], tool['validator'],
        plan['registers'], 5, len(plan['layer_offsets']) - 1, plan.get('prefix_instructions', 0),
        optimizer=tool['optimizer'], **kwargs)
    argv = [tool['disassembler'], str(destination)]; commands.append(argv)
    assembly = subprocess.check_output(argv, text=True)
    destination.with_suffix('.spvasm').write_text(assembly)
    rows[kind] = {'module': str(destination.relative_to(R)), **pin(destination),
        'original_module': pin(C / ('retained_' + stem + '.spv')),
        'registers': plan['registers'], 'program_table_equal': True, 'compiled_this_candidate': True}
    # Commit a receipt after each genuine compiler result, before metadata readers.
    write(W / 'compiler-progress.json', {'modules': rows, 'commands': commands, 'tools': tool})
    header.append('inline constexpr std::array<std::uint32_t, ' + str(len(words)) + '> k' + kind + 'NativeOwnedShader{{')
    for i in range(0, len(words), 16): header.append(','.join(str(v) + 'u' for v in words[i:i+16]) + ',')
    header.append('}};')
header.append('}')
(W / 'headers').mkdir(); (W / 'headers/retained_kernels.h').write_text('\n'.join(header) + '\n')
assert original == {str(p.relative_to(R)): pin(p) for p in C.iterdir() if p.is_file()}
write(W / 'compiler-result.json', {'modules': rows, 'commands': commands, 'tools': tool,
    'all35_original_modules_and_generated_headers_unchanged': True, 'original_inputs': original,
    'global_buffer_layout_tables_or_scalar_algorithms_changed': False, 'qualification': False})
print(json.dumps({'compiled_modules': 3, 'original35_unchanged': True, 'modules': rows}))
