"""Seal one pure-admission compile; preserve the three genuine D objects."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
basis = root / 'out/native-cpu-static-program'
tables = root / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch'
destination = root / 'bin/linux-gcc/native-cpu-dispatch-admission'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def doc(path):
    return json.loads(path.read_text())


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (work / 'prepared-batch-inputs.json').exists()
assert not (work / 'batch-compile-result.json').exists()
assert not (work / 'build/actions.json').exists()
receipts = {
    'batch-compile-result.json': 'f8434acad86353ab1a345cd35a9f4e6ef338805a2d87974ebacb49d2d9082be2',
    'build/actions.json': 'eb5d6cec32f69930ed1e426e731d7e00180d9bab53bb980d29ea9f3a121ed94c',
    'header-preparation.json': '72c76214e8b0267566801b0e7efec605c2e408f2c9ac084a892c95d7a563ef6d',
    'provider/linked-native-proof.json': '1f6259589b610272b57d4aa2fb3611a49ef5ba41c716b1b886388b7778b6d46f',
}
for name, digest in receipts.items():
    assert pin(basis / name)['sha256'] == digest, name
original = doc(basis / 'batch-compile-result.json')
assert len(original['objects']) == 3 and all(x['exit'] == 0 for x in original['objects'])
assert doc(basis / 'provider/linked-native-proof.json')['passed']
inputs = dict(original['source_inputs'])
for name, expected in inputs.items():
    assert pin(root / name) == expected, name
for item in original['objects']:
    value = item['object']
    expected = {k: value[k] for k in ('bytes', 'sha256')}
    assert pin(root / value['path']) == expected
    inputs[value['path']] = expected
guard = root / 'out/native-cpu-active-prefix/native_packet.hpp'
assert pin(guard)['sha256'] == '9df6d65d1c89016694a7f242a8fdd61979db6826c4f0d2619b2b8954a9ca6eee'
inputs[str(guard.relative_to(root))] = pin(guard)
for name in receipts:
    inputs[str((basis / name).relative_to(root))] = pin(basis / name)

source = work / 'admission.cpp'
text = source.read_text()
assert '#include "../native-cpu-active-prefix/native_packet.hpp"' in text
contracts = {}
for slug, wanted in [('transport', (230, 4699, 15967)),
                     ('endpoint', (225, 3829, 21143)),
                     ('dopri_phase', (1810, 4149, 10689))]:
    adapter = tables / (slug + '.cpp')
    table = tables / (slug + '.program.hpp')
    match = re.search(r'groups, (\d+), (\d+),\s+sirius_cpu_bridge::' + slug + r'_program, (\d+)\)',
                      adapter.read_text())
    assert match and tuple(map(int, match.groups())) == wanted
    assert re.search(r'admit_packet\(ADMISSION_VALUES, ' + str(wanted[0]) + r'U, ' + str(wanted[1]) +
                     r'U,\s+sirius_cpu_bridge::' + slug + '_program, ' + str(wanted[2]) + r'U\)', text)
    words = re.findall(r'\b\d+u\b', table.read_text())
    assert len(words) == wanted[2]
    for path in (adapter, table):
        name = str(path.relative_to(root))
        assert pin(path) == original['source_inputs'][name]
        inputs[name] = pin(path)
    contracts[slug] = dict(zip(('input_row_words', 'output_row_words', 'program_words'), wanted))
inputs[str(source.relative_to(root))] = pin(source)


def command_shape(argv):
    result = list(argv)
    result[result.index('-c') + 1] = 'SOURCE'
    result[result.index('-o') + 1] = 'OUTPUT'
    return result


template = list(original['objects'][0]['argv'])
assert all(command_shape(x['argv']) == command_shape(template) for x in original['objects'])
assert '-march=x86-64' in template
assert all(flag in template for flag in ('-std=c++2c', '-O2', '-mfpmath=sse',
    '-fexcess-precision=standard', '-fno-fast-math', '-fno-associative-math',
    '-fno-unsafe-math-optimizations', '-frounding-math', '-ffp-contract=off'))
assert str(tables) in template and str(guard.parent) in template
output = destination / 'admission.o'
assert not output.exists()
destination.mkdir(parents=True, exist_ok=True)
argv = list(template)
argv[argv.index('-c') + 1] = str(source)
argv[argv.index('-o') + 1] = str(output)
prepared = {'inputs': inputs, 'objects': [{'stage': 'PureAdmission', 'argv': argv,
    'output': str(output.relative_to(root))}], 'packet_contracts_from_original_CPP': contracts,
    'scope': 'One pure C ABI admission object only. Same pinned B guard and original canonical tables; exact D compiler/ISA/SSE/strictFP argv except source and output. Three D objects and381retainedfunctions remain unchanged; no native shader execution or qualification.'}
(work / 'prepared-batch-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
old = doc(basis / 'build/actions.json')['compile']
frozen = {name for name in old['frozen_inputs'] if name.startswith('/')}
frozen |= set(inputs)
frozen |= {str(p.relative_to(root)) for p in work.iterdir() if p.is_file()}
frozen |= {str(p.relative_to(root)) for p in (work / 'build').iterdir() if p.is_file()}
frozen |= {'out/native-cpu-slang-bridge/build_batch.py',
           'out/native-cpu-row-ownership/cold-frame/run_owned.py',
           'out/native-cpu-row-ownership/cold-frame/owned_guard.py'}
assert all((root / name).is_file() for name in frozen)
action = {**old, 'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(work / 'build_batch.py')],
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'frozen_inputs': sorted(frozen),
    'scope': 'One strict pure-admission object compile only, preserving original tools, D object bytes and source proof. No generated native-header/shader/host regeneration, provider load or native execution.'}
(work / 'build/actions.json').write_text(json.dumps({'compile': action}, indent=2) + '\n')
print(json.dumps({'prepared_objects': 1, 'preserved_D_objects': 3, 'sealed_inputs': len(frozen), 'compiled': False}))
