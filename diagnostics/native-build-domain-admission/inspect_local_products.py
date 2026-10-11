"""Read full production payloads and the actual partial-readiness observation.

Retained after the equivalent inline reader completed; this file was not
re-executed. It refuses to overwrite the completed result.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, struct, subprocess

WORK = Path(__file__).resolve().parent
HEAD = '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
header = Path('bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h')
def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'])
assert seal(header) == {'bytes': 7952166, 'sha256': '3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'}
arrays = {}
pattern = r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};'
for count, name, body in re.findall(pattern, header.read_text(), re.S):
    words = [int(value[:-1], 0) for value in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body)]
    assert len(words) == int(count) and name not in arrays
    arrays[name] = struct.pack('<' + 'I' * len(words), *words)
assert len(arrays) == 40
products = {}
for name, relative in {'app': 'bin/linux-gcc/src/sirius/app/sirius', 'backend_tests': 'bin/linux-gcc/tests/backend/sirius_backend_tests'}.items():
    path = Path(relative)
    data, before = path.read_bytes(), seal(path)
    assert data[:6] == b'\x7fELF\x02\x01'
    table = struct.unpack_from('<Q', data, 0x28)[0]
    stride, count = struct.unpack_from('<HH', data, 0x3a)
    assert stride == 64 and table + stride * count <= len(data)
    sections = []
    for index in range(count):
        fields = struct.unpack_from('<IIQQQQIIQQ', data, table + index * stride)
        kind, flags, address, offset, length = fields[1:6]
        if kind == 1 and flags & 2 and not flags & 5:
            assert offset + length <= len(data)
            sections.append({'section_index': index, 'file_offset': offset, 'address': address, 'bytes': length, 'flags': flags})
    bindings = {}
    for array, payload in arrays.items():
        matches = []
        for section in sections:
            begin = section['file_offset']
            position = data[begin:begin + section['bytes']].find(payload)
            if position >= 0:
                matches.append({**section, 'payload_file_offset': begin + position, 'payload_address': section['address'] + position})
        assert matches, (name, array)
        bindings[array] = {'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest(), 'readonly_ELF_bindings': matches}
    assert seal(path) == before
    products[name] = {'path': relative, 'executable_seal': before, 'bindings': bindings}
observation = json.loads((WORK / 'readiness-observation.json').read_text())
payload = observation['readiness_payload']
ideal = payload['ultimate_ideal']
assert payload['ready'] is False and payload['release_ready'] is False
assert payload['mandatory_build_gate']['valid'] is False
assert ideal['source_revision'] == HEAD and ideal['external_domains_admitted'] == 2
assert ideal['external_domains_required'] == 8 and len(ideal['external_domain_ids_pending']) == 6
assert ideal['satisfied'] is False
result = {'accepted': True, 'source_revision': HEAD, 'production_header': seal(header), 'production_arrays': 40,
          'readonly_complete_payload_memberships': 80, 'products': products, 'readiness_source_exact': True,
          'actual_ready': False, 'actual_release_ready': False, 'mandatory_gate_valid': False,
          'actual_admitted_domains': ideal['external_domain_ids_admitted'], 'actual_pending_domains': ideal['external_domain_ids_pending'],
          'scope': 'Whole payload/source/build and actual informational partial-receipt consumer only; no scientific dispatch/frame/full Mandatory/native runtime/release qualification.'}
destination = WORK / 'local-consumer-review.json'
assert not destination.exists()
destination.write_text(json.dumps(result, indent=2) + '\n')
