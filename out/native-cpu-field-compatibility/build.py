"""Compile only the shim, link retained/new actual objects, associate exact code."""
from pathlib import Path
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
OLD = R / 'out/native-cpu-accumulate-compatibility'
sys.path.insert(0, str(OLD))
from elf_proof_complete import pin, associate, rename_proof

assert sys.argv[1:] == ['build', str(W / 'build')]
prepared = json.loads((W / 'prepared.json').read_text())


def identity(path):
    return {'path': str(path.relative_to(R)), **pin(path)}


def unchanged():
    for name, expected in prepared['source_inputs'].items():
        assert pin(R / name) == expected, name


unchanged()
assert not (W / 'build-result.json').exists()
result = {'passed': False, 'prepared': pin(W / 'prepared.json'), 'commands': [],
          'diagnostic_shim_compilations': 0, 'new_links': 0}
for key in ('compile_argv', 'link_argv'):
    argv = prepared[key]
    assert not Path(argv[argv.index('-o')+1]).exists()
    code = subprocess.run(argv, cwd=R).returncode
    result['commands'].append({'argv': argv, 'exit': code})
    (W / 'build-result.json').write_text(json.dumps(result, indent=2) + '\n')
    if code:
        raise SystemExit(code)
    result['diagnostic_shim_compilations' if key == 'compile_argv' else 'new_links'] = 1
old = json.loads((OLD / 'compile-completion.json').read_text())
baseline = R / prepared['baseline']; executable = R / prepared['executable']
assert rename_proof(R / prepared['original_baseline_object']['path'], baseline,
    json.loads((OLD / 'prepared.json').read_text())['rename_mapping']) == json.loads((OLD / 'rename-proof.json').read_text())
associated = [identity(baseline), *prepared['objects'][:2], identity(executable.parent / 'shim.o')]
proof = associate(executable, [{**row, 'path': str(R / row['path'])} for row in associated])
proof['executable'] = identity(executable)
for row, original in zip(proof['objects'], associated):
    row['object'] = original
proof['scope'] = 'Every retained function in original renamed baseline, actual new Transport CPP/layer, diagnostic ABI shim. Only declared PC relocation operands masked; no linker relaxation. Numerical acceptance remains separate.'
(W / 'link-proof.json').write_text(json.dumps(proof, indent=2) + '\n')
resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
(W / 'build/runtime-resolution.txt').write_text(resolution)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
libraries = {str(Path(name).resolve(strict=True)) for name in re.findall(r'(/[^\s()]+)', resolution) if Path(name).is_file()}
unchanged()
result.update(passed=True, executable=identity(executable), link_proof=identity(W / 'link-proof.json'),
    driver_object=old['driver_object'], renamed_baseline_object=old['renamed_baseline_object'],
    rename_proof=old['rename_proof'], shim_object=identity(executable.parent / 'shim.o'),
    candidate_transport_object=prepared['candidate_transport_object'], source_inputs_unchanged=True,
    runtime_libraries={name: pin(Path(name)) for name in sorted(libraries)})
(W / 'build-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'built': True, 'new_shim_compilations': 1, 'new_links': 1}))
