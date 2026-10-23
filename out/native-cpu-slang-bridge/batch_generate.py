"""Generate T/E CPU source once; reuse the already authenticated DP generation."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-slang-bridge/batch-source'
authority = json.loads((R / 'out/native-cpu-row-ownership/compiler-result.json').read_text())
def pin(p):
    data = p.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
dp = R / 'bin/linux-gcc/native-cpu-slang-bridge/primitive-lowering/dopri_phase.cpp'
assert pin(dp) == json.loads((W / 'compile-result.json').read_text())['source_inputs'][str(dp.relative_to(R))]
rows = {'DopriPhase': {'source': str(dp.relative_to(R)), **pin(dp),
        'reflection': str(dp.with_suffix('.reflection.json').relative_to(R)),
        'generated_this_batch': False, 'SPIRV': authority['modules']['DopriPhase']}}
B.mkdir(parents=True, exist_ok=False)
for kind, stem in [('Transport', 'transport'), ('Endpoint', 'endpoint')]:
    old, = [a for a in authority['commands'] if a[0] == authority['tools']['compiler']
            and Path(a[1]).name == 'retained_' + stem + '.slang']
    assert (W / 'shaders' / Path(old[1]).name).read_bytes() == Path(old[1]).read_bytes()
    args = list(old)
    args[1] = str(W / 'shaders' / Path(old[1]).name)
    args[2:2] = ['-DSIRIUS_RETAINED_CPU_CPP=1', '-I', str(W / 'shaders')]
    assert args[args.index('-target') + 1] == 'spirv'
    args[args.index('-target') + 1] = 'cpp'
    index = args.index('-profile')
    assert args[index + 1] == 'spirv_1_5'
    del args[index:index + 2]
    output = B / (stem + '.cpp')
    reflection = B / (stem + '.reflection.json')
    args[args.index('-o') + 1] = str(output)
    args += ['-fp-mode', 'precise', '-reflection-json', str(reflection)]
    result = subprocess.run(args)
    row = {'original_SPIRV_compiler_argv': old, 'actual_CPP_argv': args,
           'exit': result.returncode, 'SPIRV': authority['modules'][kind],
           'generated_this_batch': True}
    if output.exists():
        row.update(source=str(output.relative_to(R)), **pin(output))
    if reflection.exists():
        row['reflection'] = str(reflection.relative_to(R))
        row['reflection_identity'] = pin(reflection)
    rows[kind] = row
    (W / 'batch-generation.json').write_text(json.dumps({'modules': rows,
        'scope': 'Source generation only; no native object or kernel execution.'}, indent=2) + '\n')
    if result.returncode:
        raise SystemExit(result.returncode)
print(json.dumps({'new_generations': 2, 'reused_DP_generation': True,
                  'source_bytes': {k: v['bytes'] for k, v in rows.items()}}))
