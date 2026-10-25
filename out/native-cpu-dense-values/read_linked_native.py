"""Retain original object/code/SPIRV joins; follow the actual eight-object graph."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-slang-bridge/read_linked_native.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '9a94c2b45471807ad74b28060edc60e19f33d42e7b9f94ba6126065d79c6f064'
source = basis.read_text()
anchors = '''    # ISPC local helper names repeat across stages. A unique exported
    # entrypoint anchors each retained object section without choosing an
    # arbitrary same-named helper from another stage.
    section_bases = {}
    for symbol in obj.symbols:
        if symbol['type'] != 2 or not symbol['size'] or not symbol['section']:
            continue
        matches = [entry for entry in linked.symbols if entry['name'] == symbol['name'] and entry['section']]
        if len(matches) == 1:
            section_bases.setdefault(symbol['section'], set()).add(matches[0]['value'] - symbol['value'])
'''
graph = '''callback = subprocess.check_output(['objdump', '-d', '--disassemble=cs_exec_fn', str(library)], text=True)
assert '<lp_native_cpu_packet_mode>' in callback
packet_mode = subprocess.check_output(['objdump', '-d', '--disassemble=lp_native_cpu_packet_mode', str(library)], text=True)
disassembly = callback + packet_mode
for slug in ('transport', 'endpoint', 'dopri_phase'):
    assert '<sirius_native_' + slug + '_admit_packet>' in packet_mode
    assert '<sirius_native_' + slug + '_rows_try_packet>' in packet_mode
    wrapper = subprocess.check_output(['objdump', '-d', '--disassemble=sirius_native_' + slug + '_rows_try_packet', str(library)], text=True)
    for target in ('sirius_native_' + slug + '_admit_packet', 'sirius_native_' + slug + '_try_packet', 'sirius_simd_' + slug + '_packet', 'sirius_native_simd_supported'):
        assert '<' + target + '>' in wrapper
    disassembly += wrapper'''
changes = [
    ("    relocations = {}", anchors + "    relocations = {}"),
    ("        final = linked.named(symbol['name'])", """        matches = [entry for entry in linked.symbols if entry['name'] == symbol['name'] and entry['section']]
        if len(matches) != 1:
            bases = section_bases.get(symbol['section'], set())
            assert len(bases) == 1, (path, symbol['section'], bases)
            expected_address = next(iter(bases)) + symbol['value']
            matches = [entry for entry in matches if entry['value'] == expected_address]
        assert len(matches) == 1, (path, symbol['name'], matches)
        final = matches[0]"""),
    ("disassembly = subprocess.check_output(['objdump', '-d', '--disassemble=cs_exec_fn', str(library)], text=True)\nfor slug in ('transport', 'endpoint', 'dopri_phase'):\n    assert '<sirius_native_' + slug + '_try_packet>' in disassembly", graph),
    ("batch = json.loads((W / 'batch-generation.json').read_text())", "batch = json.loads((R / 'out/native-cpu-slang-bridge/batch-generation.json').read_text())"),
    ("(P / 'linked-native.disassembly.txt').write_text(all_native_disassembly)",
     "# Full native disassembly is checked in memory; avoid retaining a duplicate large dump."),
    ("'callback_direct_calls_all3': True", "'callback_shared_mode_reaches_rows_pure_SIMD_and_scalar_all3': True"),
    ('Callback includes direct native calls and ordinary fallback under existing IEEE32 guard.',
     'Callback shared guard reaches all3 rows wrappers, each retaining pure admission, SIMD8 and original scalar tails; ordinary fallback and IEEE32 scope retained. Same-named ISPC locals are associated by unique exported object-section anchors.'),
]
for before, after in changes:
    assert source.count(before) == 1, before
    source = source.replace(before, after)
exec(compile(source, str(basis), 'exec'))

# The source/object joins do not alone prove relocated call destinations.
import json
import re
import subprocess
R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
receipt = json.loads((P / 'provider-build-result.json').read_text())
library = R / receipt['library']
assembly = subprocess.check_output(['/usr/bin/objdump', '-d', str(library)], text=True)
targets = {}
for prefix, target in (('RTAddStore_0___', 'sirius_native_dense_add8'),
                       ('RTMultiplyStore_0___', 'sirius_native_dense_multiply8')):
    callers = []
    for block in re.split(r'\n(?=[0-9a-f]+ <)', assembly):
        match = re.match(r'([0-9a-f]+) <([^>]+)>:', block)
        if match and match[2].startswith(prefix):
            assert '<' + target + '>' in block, (match[2], target)
            callers.append({'symbol': match[2], 'address': int(match[1], 16), 'actual_target': target})
    assert len(callers) == 3, (prefix, callers)
    targets[target] = callers
proof = json.loads((P / 'linked-native-proof.json').read_text())
proof['dense_C_ABI_actual_targets_all3'] = targets
(P / 'linked-native-proof.json').write_text(json.dumps(proof, indent=2) + '\n')
print(json.dumps({'actual_dense_targets_joined': 6, 'stages': 3}))

# Verify new value ABI calls within each independently associated raw object.
blocks = {}
for block in re.split(r'\n(?=[0-9a-f]+ <)', assembly):
    match = re.match(r'([0-9a-f]+) <([^>]+)>:', block)
    if match:
        blocks[int(match[1], 16)] = block
value_calls = {}
for module in proof['objects'][4:7]:
    callers = {}
    for target in ('sirius_native_dense_values_add8', 'sirius_native_dense_values_multiply8'):
        found = [row for row in module['functions'] if '<' + target + '>' in blocks.get(row['linked_address'], '')]
        assert found, (module['object'], target)
        callers[target] = [{'symbol':row['symbol'], 'address':row['linked_address']} for row in found]
    value_calls[module['object']['path']] = callers
for wrapper, target in (('sirius_native_dense_values_add8','sirius_native_dense_add8'),
                        ('sirius_native_dense_values_multiply8','sirius_native_dense_multiply8')):
    row, = [row for row in proof['objects'][-1]['functions'] if row['symbol']==wrapper]
    assert '<' + target + '>' in blocks[row['linked_address']]
proof['dense_value_ABI_actual_targets_all3'] = value_calls
proof['packing_wrappers_call_unchanged_AC_math'] = True
(P / 'linked-native-proof.json').write_text(json.dumps(proof, indent=2) + '\n')
print(json.dumps({'value_ABI_target_graph_all3': True, 'unchanged_AC_math_targets': 2}))
