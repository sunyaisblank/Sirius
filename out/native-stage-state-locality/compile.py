from pathlib import Path
import hashlib, importlib.util, json, sys

sys.dont_write_bytecode = True
R = Path.cwd()
W = Path(__file__).resolve().parent
def module(name, p):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m
b = module('kernel_builder', R / 'scripts/build-retained-kernels.py')
m = module('kernel_program', R / 'src/sirius/kernels/retained_program.py')
general = m.build_transport_program(parallel=True)
special = m.build_schwarzschild_transport_program(parallel=True)
def encode(p): return [p['instructions'], p['registers'], 40] + p['outputs'] + p['operations'] + p['layer_offsets']
tables = encode(general) + encode(special)
frame = (len(encode(general)), len(encode(special)), len(special['layer_offsets']) - 1)
assert general['registers'] == special['registers'] == 459
assert (459 * 5 + 25 * 5 + 40 * 5) * 4 + 8 == 10488 <= 16384
tools = R / 'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
records = {}
for label, wide in [('candidate-default', False), ('excluded-wide', True)]:
    destination = W / 'compiler' / (label + '.spv')
    b.compile_shader(W / 'source/retained_transport.slang', destination,
                     '/opt/slang/bin/slangc', str(tools / 'spirv-as'),
                     str(tools / 'spirv-dis'), str(tools / 'spirv-val'),
                     general['registers'], 5, len(general['layer_offsets']) - 1,
                     fp64=wide, transport_specialized=frame)
    data = destination.read_bytes()
    records[label] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    if wide:
        assert data == (R / 'bin/linux-gcc/src/sirius/backend/retained/retained_transport_fp64.spv').read_bytes()
records['program_words'] = len(tables)
records['program_sha256'] = hashlib.sha256(b''.join(x.to_bytes(4, 'little') for x in tables)).hexdigest()
records['source_only_shared_bound_bytes'] = 10488
records['arithmetic_or_benefit_admission'] = False
(W / 'compiler-review.json').write_text(json.dumps(records, indent=2) + '\n')
print(json.dumps(records))
