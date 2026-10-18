import hashlib
import json
import struct
import sys
from pathlib import Path


def messagepack(data):
    offset = 0

    def take(count):
        nonlocal offset
        assert 0 <= count <= len(data) - offset
        part = data[offset:offset + count]
        offset += count
        return part

    def integer(count, signed=False):
        return int.from_bytes(take(count), 'big', signed=signed)

    def read():
        tag = integer(1)
        if tag < 0x80:
            return tag
        if tag >= 0xe0:
            return tag - 256
        if 0x80 <= tag < 0x90:
            return mapping(tag & 15)
        if 0x90 <= tag < 0xa0:
            return [read() for _ in range(tag & 15)]
        if 0xa0 <= tag < 0xc0:
            return take(tag & 31).decode('utf-8')
        if tag in (0xc0, 0xc2, 0xc3):
            return {0xc0: None, 0xc2: False, 0xc3: True}[tag]
        if 0xcc <= tag <= 0xcf:
            return integer(1 << (tag - 0xcc))
        if 0xd0 <= tag <= 0xd3:
            return integer(1 << (tag - 0xd0), True)
        if 0xd9 <= tag <= 0xdb:
            return take(integer(1 << (tag - 0xd9))).decode('utf-8')
        if tag in (0xdc, 0xdd):
            return [read() for _ in range(integer(2 if tag == 0xdc else 4))]
        if tag in (0xde, 0xdf):
            return mapping(integer(2 if tag == 0xde else 4))
        raise ValueError('Unsupported MessagePack tag %02x' % tag)

    def mapping(count):
        result = {}
        for _ in range(count):
            key, value = read(), read()
            assert key not in result
            result[key] = value
        return result

    value = read()
    assert offset == len(data)
    return value


folder = Path(sys.argv[1])
path = folder / 'shader/shader-info-1.bin'
data = path.read_bytes()
assert data[:6] == b'\x7fELF\x02\x01'
header = struct.unpack_from('<16sHHIQQQIHHHHHH', data)
assert header[2] == 224 and header[11] == 64
section_offset, section_count, names_index = header[6], header[12], header[13]
headers = [struct.unpack_from('<IIQQQQIIQQ', data, section_offset + i * 64)
           for i in range(section_count)]
names_header = headers[names_index]
names = data[names_header[4]:names_header[4] + names_header[5]]
sections = {}
metadata = None
for section in headers:
    name = names[section[0]:].split(b'\0', 1)[0].decode()
    sections[name] = {'offset': section[4], 'bytes': section[5], 'flags': section[2]}
    if section[1] != 7:
        continue
    part = data[section[4]:section[4] + section[5]]
    offset = 0
    while offset < len(part):
        name_bytes, description_bytes, kind = struct.unpack_from('<III', part, offset)
        offset += 12
        owner = part[offset:offset + name_bytes].rstrip(b'\0')
        offset += (name_bytes + 3) & ~3
        description = part[offset:offset + description_bytes]
        offset += (description_bytes + 3) & ~3
        if owner == b'AMDGPU' and kind == 32:
            assert metadata is None
            metadata = messagepack(description)
    assert offset == len(part)
assert metadata is not None
isa = (folder / 'shader/shader-info-2.isa.txt').read_bytes()
result = {'binary_bytes': len(data), 'binary_sha256': hashlib.sha256(data).hexdigest(),
          'sections': sections, 'PAL_metadata': metadata,
          'disassembly_bytes': len(isa), 'disassembly_sha256': hashlib.sha256(isa).hexdigest(),
          'scope': 'Metadata of this diagnostic pipeline only; no dispatch, dynamic occupancy, bottleneck attribution or performance claim'}
(folder / 'decoded-driver-object.json').write_text(json.dumps(result, indent=2) + '\n')
stage = metadata['amdpal.pipelines'][0]['.hardware_stages']['.cs']
print(json.dumps({'folder': str(folder), 'text_bytes': sections['.text']['bytes'], 'stage': stage}, indent=2))
