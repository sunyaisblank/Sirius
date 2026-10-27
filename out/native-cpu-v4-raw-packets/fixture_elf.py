"""Narrow original opcode association extension for two fixture GOT fields."""
from pathlib import Path
import hashlib
import struct
from elf_proof import Elf, pin

R = Path.cwd()
basis = R / 'out/native-cpu-accumulate-compatibility/elf_proof.py'
assert pin(basis)['sha256'] == '278766535b0ddcb9ad2ab97b68805cd851dda31036567a35f8bf5a07afedb0f7'
before = 'assert kind in (2, 4, 42) and offset + 4 <= symbol[\'size\']'
after = 'assert kind in (2, 4, 9, 42) and offset + 4 <= symbol[\'size\']'
assert basis.read_text().count(before) == 1
field = R / 'out/native-cpu-field-compatibility/elf_proof.py'
assert pin(field)['sha256'] == 'faf7dd9d657941df283a36cf9f3b553a066ef18e4e46626689804c35e0ec6cd7'
source = field.read_text()
assert source.count('source = basis.read_text()') == 1
source = source.replace('source = basis.read_text()',
    'source = basis.read_text().replace(' + repr(before) + ', ' + repr(after) + ', 1)', 1)
# Keep its existing unique-local-section anchors and dynamic-table support.
scope = {}
exec(compile(source, str(field), 'exec'), scope)


def word(elf, address, fmt):
    size = struct.calcsize(fmt)
    sections = [(i, section) for i, section in enumerate(elf.sections)
                if section[2] & 2 and section[1] != 8 and
                section[3] <= address and address + size <= section[3] + section[5]]
    i, section = sections[0]
    assert len(sections) == 1
    return struct.unpack_from(fmt, elf.bytes(i), address - section[3])[0]


def associate(executable, objects):
    linked = Elf(executable)
    joins = []
    for item in objects:
        obj = Elf(Path(item['path']))
        for section, offset, kind, addend, target in obj.relocations:
            if kind != 9:
                continue
            # This extension is for the unchanged Transport fixture's two
            # indirect raw calls. Every other use needs a separate review.
            assert Path(item['path']).name == 'transport-driver.o'
            assert target['name'] in {'baseline_sirius_simd_transport_packet',
                                       'sirius_simd_v4_transport_packet'}
            assert addend == -4
            main, = [s for s in obj.symbols if s['name'] == 'main' and s['section']]
            assert main['section'] == section and main['value'] <= offset < main['value'] + main['size']
            final_main, = [s for s in linked.symbols if s['name'] == 'main' and s['section']]
            field = final_main['value'] + offset - main['value']
            slot = field + word(linked, field, '<i') - addend
            destination, = [s for s in linked.symbols if s['name'] == target['name'] and s['section']]
            actual = word(linked, slot, '<Q')
            assert actual == destination['value']
            dynamic, = [(kind, addend, target) for _, location, kind, addend, target
                         in linked.relocations if location == slot]
            assert dynamic[0] == 8 and dynamic[1] == destination['value'] and not dynamic[2]['name']
            joins.append(dict(object=item['path'], relocation_type=9,
                target=target['name'], field_address=field, GOT_slot=slot,
                actual_address=actual, expected_address=destination['value'],
                runtime_relative_relocation_addend_checked=True))
    assert len(joins) in (0, 2)
    result = scope['associate'](executable, objects)
    result['declared_GOTPCREL_fields_resolve_exact_targets'] = joins
    result['scope'] += ' Two original Transport main GOTPCREL fields additionally resolve to exact original baseline/new-v4 functions; no other type9 admitted.'
    return result
