"""Generate one original DP shader as C++ source; never compile or execute it."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-slang-bridge'
authority = json.loads((R / 'out/native-cpu-row-ownership/compiler-result.json').read_text())
original, = [a for a in authority['commands'] if a[0] == authority['tools']['compiler']
             and Path(a[1]).name == 'retained_dopri_phase.slang']
args = list(original)
assert args[args.index('-target') + 1] == 'spirv'
args[args.index('-target') + 1] = 'cpp'
index = args.index('-profile')
assert args[index + 1] == 'spirv_1_5'
del args[index:index + 2]
B.mkdir(parents=True, exist_ok=True)
output = B / 'dopri_phase.cpp'
reflection = B / 'dopri_phase.reflection.json'
assert not output.exists() and not reflection.exists()
args[args.index('-o') + 1] = str(output)
args += ['-fp-mode', 'precise', '-reflection-json', str(reflection)]
result = subprocess.run(args)
def pin(p):
    data = p.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
record = {'original_SPIRV_compiler_argv': original, 'actual_CPP_argv': args,
          'exit': result.returncode, 'scope': 'One C++ source generation only; no native object, provider or numerical/performance execution.'}
for p in [output, reflection]:
    if p.exists():
        record[str(p.relative_to(R))] = pin(p)
(W / 'generation-result.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
raise SystemExit(result.returncode)
