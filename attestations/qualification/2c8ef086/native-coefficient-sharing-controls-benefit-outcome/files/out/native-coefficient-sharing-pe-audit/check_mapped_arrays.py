"""Separate ROOT implementation: follow MAP addresses into each PE directly.

The mature audit independently searches full payloads before joining MAP RVAs.
This reader starts from MAP symbols and never imports that implementation.
"""
import hashlib
import json
import re
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
STAGE = ROOT / 'bin/windows-msvc/native-coefficient-sharing-matched-current'
INPUTS = ROOT / 'out/native-coefficient-sharing-producer-inputs'


def signature(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def header_arrays(path):
    source = path.read_text()
    arrays = {}
    for count, name, data in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', source, re.S):
        words = [int(value) for value in re.findall(r'(\d+)u', data)]
        assert len(words) == int(count) and name not in arrays
        arrays[name] = struct.pack('<%dI' % len(words), *words)
    return arrays


record = json.loads((STAGE / 'producer-result.json').read_text())
assert record['source_revision'] == 'e1c618f83861c478fe18c5a22f5d9aa27f6280fd'
assert record['passed'] and record['run_id'] == '39393948651' and record['run_attempt'] == '1'
results = {}
decoded = {}
for arm in ('baseline', 'candidate'):
    path = STAGE / 'payload' / arm / ('original_' + arm + '_controls.exe')
    data = path.read_bytes()
    assert data[:2] == b'MZ'
    pe = struct.unpack_from('<I', data, 60)[0]
    assert data[pe:pe + 4] == b'PE\0\0'
    machine, count = struct.unpack_from('<HH', data, pe + 4)
    optional_size = struct.unpack_from('<H', data, pe + 20)[0]
    optional = pe + 24
    assert machine == 0x8664 and struct.unpack_from('<H', data, optional)[0] == 0x20b
    image_base = struct.unpack_from('<Q', data, optional + 24)[0]
    sections = []
    for index in range(count):
        offset = optional + optional_size + 40 * index
        name, virtual_size, rva, raw_size, file_offset = struct.unpack_from('<8sIIII', data, offset)
        assert file_offset + raw_size <= len(data)
        sections.append((rva, raw_size, file_offset, name.rstrip(b'\0').decode()))
    mapping = path.with_suffix('.map').read_text()
    mapped = {}
    for line in mapping.splitlines():
        match = re.match(r'\s*[0-9a-fA-F]{4}:[0-9a-fA-F]+\s+(\S+)\s+([0-9a-fA-F]{8,16})\b', line)
        if match:
            mapped.setdefault(match[1], []).append(int(match[2], 16) - image_base)
    arm_result = {}
    decoded[arm] = {}
    for filename, namespace, expected_count in [('retained_kernels.h', 'retained_program', 40),
                                                 ('retained_owner_kernels.h', 'retained_owner', 2)]:
        arrays = header_arrays(INPUTS / arm / filename)
        assert len(arrays) == expected_count
        decoded[arm][namespace] = arrays
        for name, expected in arrays.items():
            prefix = '?' + name + '@' + namespace + '@backend@sirius@@'
            addresses = {rva for symbol, values in mapped.items() if symbol.startswith(prefix) for rva in values}
            assert addresses, (arm, name)
            joins = []
            for address in sorted(addresses):
                selected = [(rva, raw_size, file_offset, section) for rva, raw_size, file_offset, section in sections
                            if rva <= address and address + len(expected) <= rva + raw_size]
                assert len(selected) == 1
                rva, raw_size, file_offset, section = selected[0]
                begin = file_offset + address - rva
                assert data[begin:begin + len(expected)] == expected, (arm, name, address)
                joins.append({'rva': address, 'file_offset': begin, 'section': section})
            arm_result[namespace + '::' + name] = {**signature(expected), 'MAP_to_PE_joins': joins}
    assert len(arm_result) == 42
    assert sum(symbol.startswith('?Create@RetainedCompute@') for symbol in mapped) == 1
    assert not re.search(r'sirius_backend(?:\.lib)?[:(][^\s)]*retained_compute\.obj', mapping, re.I)
    results[arm] = {'PE': signature(data), 'arrays': arm_result}
main_old, main_new = decoded['baseline']['retained_program'], decoded['candidate']['retained_program']
assert {name for name in main_old if main_old[name] != main_new[name]} == {'kTransportShader'}
assert decoded['baseline']['retained_owner'] == decoded['candidate']['retained_owner']
report = {'schema': 'sirius-native-coefficient-sharing-separate-ROOT-producer-review-v1',
          'passed': True, 'producer_source_revision': record['source_revision'],
          'run_id': record['run_id'], 'reviewer': 'ROOT using a separate stdlib implementation; no independent-agent claim',
          'complete_mapped_array_checks': 84, 'main_changed': ['kTransportShader'], 'main_exact': 39,
          'owner2_exact': True, 'one_private_RetainedCompute_Create_per_arm': True,
          'results': results, 'runtime_admitted': False,
          'limits': 'Actual compiled payload/MAP/address check only. Mature source/build/tool/CL/LINK/canonical archive audit remains separately required. No numeric, benefit, science, frame or full qualification claim.'}
output = WORK / 'separate-ROOT-producer-review.json'
assert not output.exists()
output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: value for key, value in report.items() if key != 'results'}))
