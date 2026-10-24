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
    for target in ('sirius_native_' + slug + '_admit_packet', 'sirius_native_' + slug + '_try_packet', 'sirius_simd_' + slug + '_packet'):
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
