import json
import re
import struct

def arrays(path):
    result = {}
    for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', path.read_text(), re.S):
        words = [int(v) for v in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count) and name not in result
        result[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    return result

def elf_symbols(path):
    data = path.read_bytes()
    h = struct.unpack_from('<16sHHIQQQIHHHHHH', data)
    assert h[0][:6] == b'\x7fELF\x02\x01' and h[2] == 62
    assert h[11] == 64
    sections = [struct.unpack_from('<IIQQQQIIQQ', data, h[6] + i*h[11]) for i in range(h[12])]
    result = {}
    for section in sections:
        if section[1] != 2:
            continue
        strings = sections[section[6]]
        strings_data = data[strings[4]:strings[4]+strings[5]]
        assert section[9] == 24
        for off in range(section[4], section[4]+section[5], section[9]):
            name, info, other, index, value, size = struct.unpack_from('<IBBHQQ', data, off)
            label = strings_data[name:strings_data.find(b'\0', name)].decode()
            if label and 0 < index < len(sections):
                sec = sections[index]
                if label.startswith(('_ZN6sirius7backend16retained_program', '_ZN6sirius7backend15RetainedCompute6Create', '_ZN6sirius7backend14retained_owner')):
                    assert label not in result, ('duplicate required ELF symbol', label)
                raw = None
                if sec[1] != 8 and value >= sec[3] and value + size <= sec[3]+sec[5]:
                    raw = data[sec[4]+value-sec[3]:sec[4]+value-sec[3]+size]
                result[label] = {'address': value, 'size': size, 'binding': info >> 4,
                                 'type': info & 15, 'section_flags': sec[2], 'data': raw}
    return result

def elf_arrays(path, expected, namespace='retained_program'):
    symbols = elf_symbols(path)
    prefix = '_ZN6sirius7backend' + str(len(namespace)) + namespace
    found = {name: symbols[prefix + str(len(name)) + name + 'E'] for name in expected}
    records = {}
    for name, symbol in found.items():
        assert symbol['type'] == 1 and symbol['section_flags'] & 2, name
        assert symbol['size'] == len(expected[name]) and symbol['data'] == expected[name], name
        records[name] = {'bytes': symbol['size'], 'sha256': hashlib.sha256(symbol['data']).hexdigest(),
                         'address': symbol['address'], 'binding': symbol['binding']}
    return records, symbols
