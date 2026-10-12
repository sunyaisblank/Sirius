"""ELF section, function, symbol and symbolic relocation reader."""
import struct

def elf(path, label_prefix=".LBB"):
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
            if name.startswith(label_prefix):
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
