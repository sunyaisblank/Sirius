"""Require complete code, constants and symbolic relocation correspondence."""
from pathlib import Path
import hashlib
import json
import struct

root = Path.cwd()
work = Path(__file__).resolve().parent
folder = root / 'bin/linux-gcc/portable-execution-attribution'
reference = root / 'attestations/software-vulkan/b5d2abf/jit-stage-linked-byte-binding/capture/jit-10.elf'
sha = lambda data: hashlib.sha256(data).hexdigest()

def elf(path):
    data = path.read_bytes()
    assert data[:6] == b'\x7fELF\x02\x01'
    assert struct.unpack_from('<HH', data, 16) == (1, 62)
    offset = struct.unpack_from('<Q', data, 40)[0]
    width, count, names_index = struct.unpack_from('<HHH', data, 58)
    assert width == 64 and offset + count * width <= len(data)
    sections = []
    for index in range(count):
        values = struct.unpack_from('<IIQQQQIIQQ', data, offset + index * width)
        sections.append(dict(zip(('name_at','type','flags','address','offset','bytes',
                                  'link','info','align','entry_bytes'), values)))
    names = sections[names_index]
    strings = data[names['offset']:names['offset'] + names['bytes']]
    for section in sections:
        at = section['name_at']
        section['name'] = strings[at:strings.index(b'\0', at)].decode()
        if section['type'] != 8:
            assert section['offset'] + section['bytes'] <= len(data)
        section['data'] = data[section['offset']:section['offset'] + section['bytes']]
    tables, functions, labels = {}, {}, {}
    for index, section in enumerate(sections):
        if section['type'] != 2:
            continue
        assert section['entry_bytes'] == 24
        strings = sections[section['link']]['data']
        table = []
        for at in range(0, section['bytes'], 24):
            name, info, other, target, value, size = struct.unpack_from('<IBBHQQ', section['data'], at)
            name = strings[name:strings.index(b'\0', name)].decode()
            symbol = {'name': name, 'type': info & 15, 'section': target, 'value': value, 'size': size}
            table.append(symbol)
            if target and target < count and info & 15 == 2:
                assert name not in functions
                assert value + size <= sections[target]['bytes']
                functions[name] = (sections[target]['name'], value, size)
            if name.startswith('.LBB'):
                assert name not in labels and 0 < target < count
                labels[name] = (sections[target]['name'], value)
        tables[index] = table
    relocations = {}
    for section in sections:
        if section['type'] != 4:
            continue
        assert section['entry_bytes'] == 24
        target_name = sections[section['info']]['name']
        records = []
        for at in range(0, section['bytes'], 24):
            where, info, addend = struct.unpack_from('<QQq', section['data'], at)
            symbol = tables[section['link']][info >> 32]
            kind = info & 0xffffffff
            if symbol['section'] == 0:
                target = ('undefined', symbol['name'], addend)
            elif symbol['section'] == 0xfff1:
                target = ('absolute', symbol['name'], symbol['value'] + addend)
            else:
                assert symbol['section'] < count
                # Local/section symbol spellings can differ with label retention;
                # their complete defined-target section-relative expression must agree.
                target = ('defined', sections[symbol['section']]['name'], symbol['value'] + addend)
            records.append((where, kind, target))
        assert target_name not in relocations
        relocations[target_name] = sorted(records)
    allocated = {s['name']: s for s in sections if s['flags'] & 2 and s['bytes']}
    return data, allocated, functions, relocations, labels

before = elf(reference)
assert set(before[2]) == {'cs_variant', 'cs_co_variant', 'cs_co_variant.resume', 'cs_co_variant.destroy'}
reports = []
for name in ('endpoint.o', 'endpoint-llvm-labels.o'):
    current = elf(folder / name)
    assert before[2] == current[2], (name, 'whole function layout')
    assert set(before[1]) == set(current[1]), (name, 'allocated section set')
    joins = []
    for section, original in before[1].items():
        after = current[1][section]
        assert (original['flags'], original['align'], original['bytes']) == (
            after['flags'], after['align'], after['bytes']), (name, section, 'layout')
        assert original['data'] == after['data'], (name, section, 'complete bytes')
        assert before[3].get(section, []) == current[3].get(section, []), (name, section, 'symbolic relocation records')
        joins.append({'section': section, 'bytes': original['bytes'], 'sha256': sha(original['data']),
                      'symbolic_relocations': len(before[3].get(section, []))})
    reports.append({'artifact': name, 'artifact_sha256': sha(current[0]),
                    'all_four_complete_functions_layout_and_bytes_exact': True,
                    'all_allocated_sections_and_symbolic_relocations_exact': joins})
labels = elf(folder / 'endpoint-llvm-labels.o')[4]
assert labels and all(section == '.ltext' and 0 <= value < before[1]['.ltext']['bytes']
                      for section, value in labels.values())
(work / 'basic-block-label-offsets.json').write_text(json.dumps(labels, indent=2) + '\n')
result = {'passed': True, 'reference_sha256': sha(before[0]), 'reports': reports,
          'complete_functions_each': 4, 'complete_allocated_sections_each': len(before[1]),
          'retained_machine_basic_block_labels': len(labels),
          'producer_options_uniquely_recovered': False,
          'scope': 'Exact offline artifact offset correspondence only. IR block comments '
                   'may reflect merging/splitting/motion; no universal original-source map, '
                   'current live PCs, exclusive costs, speed or fullqualification claim.'}
(work / 'emission-verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
