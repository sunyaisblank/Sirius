import hashlib
import json
from pathlib import Path
import re
import struct

root = Path.cwd()
work = Path(__file__).resolve().parent
build = json.loads((work/'restoration-build-and-arrays.json').read_text())
header = root/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert hashlib.sha256(header.read_bytes()).hexdigest() == build['whole_header_sha256']
payloads = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', header.read_text(), re.S):
    words = [int(q[:-1], 0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert len(words) == int(count)
    payloads[name] = struct.pack('<'+'I'*len(words), *words)
assert len(payloads) == 40
consumers = {}
for relative in ['bin/linux-gcc/src/sirius/app/sirius', 'bin/linux-gcc/tests/backend/sirius_backend_tests', 'bin/linux-gcc/src/sirius/app/sirius_app_tests', 'bin/linux-gcc/src/sirius/app/sirius_render_tests']:
    data = (root/relative).read_bytes()
    assert data[:6] == b'\x7fELF\x02\x01'
    offset = struct.unpack_from('<Q', data, 40)[0]
    size, count, name_index = struct.unpack_from('<HHH', data, 58)
    assert size == 64 and 0 < count < 10000
    sections = [struct.unpack_from('<IIQQQQIIQQ', data, offset+i*size) for i in range(count)]
    names = sections[name_index]
    strings = data[names[4]:names[4]+names[5]]
    bindings = {}
    for name, payload in payloads.items():
        found = []
        for section in sections:
            # PROGBITS, allocated, without writable or executable flags.
            if section[1] != 1 or not section[2] & 2 or section[2] & 5:
                continue
            body = data[section[4]:section[4]+section[5]]
            at = body.find(payload)
            if at >= 0:
                found.append(dict(section=strings[section[0]:].split(b'\0', 1)[0].decode(),
                                  flags=section[2], file_offset=section[4]+at))
        assert found, (relative, name)
        identity = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
        assert identity == build['arrays'][name]
        bindings[name] = dict(**identity, occurrences=found)
    consumers[relative] = dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), arrays=bindings)
record = dict(source_revision=build['source_revision'], whole_header_sha256=build['whole_header_sha256'],
              actual_consumers=consumers, pass_=True,
              scope='All40 complete payloads in allocated nonwritable/nonexecutable ELF sections of all four actual Linux consumers. No all40 loaded/runtime-selection or native/performance/fullqualification claim.')
(work/'restoration-linux-readonly-payloads.json').write_text(json.dumps(record, indent=2)+'\n')
print('All40 complete actual Linux consumer payloads bound in allocated read-only sections.')
