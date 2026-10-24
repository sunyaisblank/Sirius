"""Seal one medium-tail adapter compile with the exact genuine I C++ argv."""
from pathlib import Path
import hashlib
import json
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
original = root / 'out/native-cpu-simd-packet-policy'
previous = root / 'out/native-cpu-nir-dead-progress'


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (work / 'prepared-adapter.json').exists()
assert not (work / 'adapter-compile-result.json').exists()
assert not (work / 'compile/actions.json').exists()
compiled = json.loads((original / 'adapter-compile-result.json').read_text())
owner = json.loads((original / 'compile/adapter/owner.json').read_text())
assert compiled['all_prepared_source_bytes_unchanged'] and owner['passed'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes']
item, = compiled['objects']
assert item['exit'] == 0
argv = list(item['argv'])
assert argv[0] == '/usr/bin/x86_64-linux-gnu-g++-14'
assert argv[argv.index('-c') + 1] == str(original / 'rows_adapter.cpp')
assert all(flag in argv for flag in ('-std=c++2c', '-O2', '-march=x86-64', '-mfpmath=sse',
    '-fno-fast-math', '-fno-associative-math', '-fno-unsafe-math-optimizations',
    '-frounding-math', '-ffp-contract=off'))
contract = json.loads((previous / 'provider/source-contract.json').read_text())
assert contract['native_objects'][7] == item['object']
proof = json.loads((work / 'partition-proof.json').read_text())
assert proof['adapter_only_three_group_predicates_changed']
assert proof['seven_reused_native_objects'] == contract['native_objects'][:7]
inputs = {**compiled['source_inputs'], **contract['source_inputs']}
# I's result records its actual argv, while D's original owned compile
# snapshot supplies the original compiler/runtime/header byte identities.
tool_snapshot = root / 'out/native-cpu-static-program/build/compile/inputs-before.json'
tools = json.loads(tool_snapshot.read_text())
tool_key = str(Path(argv[0]).resolve(strict=True))
assert tool_key in tools
for name, value in tools.items():
    path = Path(name)
    if path.is_absolute() and not path.is_relative_to(root):
        inputs[name] = {key: value[key] for key in ('bytes', 'sha256')}
inputs[str(tool_snapshot.relative_to(root))] = pin(tool_snapshot)
for value in contract['native_objects']:
    inputs[value['path']] = {key: value[key] for key in ('bytes', 'sha256')}
for path in (original / 'adapter-compile-result.json', original / 'compile/adapter/owner.json',
             previous / 'provider/source-contract.json'):
    inputs[str(path.relative_to(root))] = pin(path)
for path in work.rglob('*'):
    if path.is_file():
        inputs[str(path.relative_to(root))] = pin(path)
assert inputs[tool_key] == {key: tools[tool_key][key] for key in ('bytes', 'sha256')}
for name, value in inputs.items():
    assert pin(root / name) == value, name
output = root / 'bin/linux-gcc/native-cpu-simd-medium-tails/objects/rows_adapter.o'
assert not output.exists()
output.parent.mkdir(parents=True, exist_ok=True)
argv[argv.index('-c') + 1] = str(work / 'rows_adapter.cpp')
argv[argv.index('-o') + 1] = str(output)
inverse = list(argv)
inverse[inverse.index('-c') + 1] = item['argv'][item['argv'].index('-c') + 1]
inverse[inverse.index('-o') + 1] = item['argv'][item['argv'].index('-o') + 1]
assert inverse == item['argv']
prepared = {'inputs': inputs, 'argv': argv, 'output': str(output.relative_to(root))}
(work / 'prepared-adapter.json').write_text(json.dumps(prepared, indent=2) + '\n')
base = json.loads((original / 'compile/actions.json').read_text())['adapter']
assert base['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
action = {**base, 'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(work / 'build_adapter.py')],
    'frozen_inputs': sorted(set(inputs) | {str((work / 'prepared-adapter.json').relative_to(root))}),
    'scope': 'N one C ABI adapter compile only, exact genuine I g++14 baseline ISA/strict argv exceptsource/output paths. Only3 group predicates admit masked4..8rows,1..3 scalar tails retained. All7 D/F/L objects unchanged; no ISPC/native math recompilation/provider/shader/workload execution.'}
(work / 'compile/actions.json').write_text(json.dumps({'adapter': action}, indent=2) + '\n')
print(json.dumps({'prepared_adapter_objects': 1, 'unchanged_objects': 7, 'compiled': False}))
