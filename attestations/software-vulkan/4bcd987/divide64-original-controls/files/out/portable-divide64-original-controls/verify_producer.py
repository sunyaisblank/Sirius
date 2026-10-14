"""Bind all forty generated arrays to actual ELF bytes and fresh compile inputs."""
from pathlib import Path
import hashlib, json, re, shlex, struct, subprocess

R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc'
T = B / 'portable-divide64-original-controls/controls'
def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def arrays(path):
    result = {}
    pattern = r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};'
    for count, name, body in re.findall(pattern, path.read_text(), re.S):
        words = [int(s) for s in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count)
        result[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    assert len(result) == 40
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
actions = json.loads((W / 'actions.json').read_text())
compile_records = {}
for label, name in [('compile-compute-sealed', 'compute-sealed'), ('compile-observer', 'observer')]:
    owner = json.loads((W / label / 'owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and owner['remaining_owned_processes'] == []
    before = json.loads((W / label / 'inputs-before.json').read_text())
    assert before == json.loads((W / label / 'inputs-after.json').read_text())
    dependencies = shlex.split((T / (name + '.d')).read_text().replace('\\\n', '').split(':', 1)[1])
    normalized = {str((R / key).resolve()): pin for key, pin in before.items()}
    records = {}
    for dependency in dependencies:
        path = (R / dependency).resolve()
        assert str(path) in normalized, ('actual dependency missing precompile seal', label, str(path))
        pin = {key: normalized[str(path)][key] for key in ['bytes', 'sha256']}
        assert identity(path) == pin, dependency
        records[str(path.relative_to(R))] = pin
    assert 'bin/linux-gcc/portable-divide64-original-controls/controls/headers/retained_kernels.h' in records
    if label == 'compile-compute-sealed':
        for header in ['retained_owner_kernels.h', 'retained_owner_layout.h']:
            assert 'bin/linux-gcc/portable-divide64-original-controls/controls/headers/' + header in records
    assert 'src/sirius/backend/vulkan/vulkan_device.h' not in records or records['src/sirius/backend/vulkan/vulkan_device.h']['sha256'] == '604c1ca5ca914d73e0e1c7288f832309e39a07b1ac3006b98812d0fa7e80f3de'
    compile_records[label] = {'owner': identity(W / label / 'owner.json'), 'object': identity(T / (name + '.o')),
                              'actual_dependencies': records}
link = json.loads((W / 'link-candidate/owner.json').read_text())
assert link['passed'] and link['status'] == 'terminal' and link['child_exit'] == 0
assert link['whole_inputs_unchanged'] and link['source_unchanged']
linkmap = (T / 'candidate.map').read_text()
assert 'libsirius_backend.a(retained_compute.cpp.o)' not in linkmap
assert str(T / 'compute-sealed.o') in linkmap and str(T / 'observer.o') in linkmap
assert '/audited-baseline-observer.o' not in linkmap
stock = arrays(B / 'src/sirius/backend/retained/retained_kernels.h')
candidate = arrays(T / 'headers/retained_kernels.h')
array_records = {}
selector_records = {}
for mode, expected in [('baseline', stock), ('candidate', candidate)]:
    executable = T / ('original_' + mode + '_controls')
    array_records[mode], symbols = elf_arrays(executable, expected)
    create_symbols = [name for name in symbols if name.startswith('_ZN6sirius7backend15RetainedCompute6Create') and '.' not in name]
    assert len(create_symbols) == 1
    create = symbols[create_symbols[0]]
    assembly = subprocess.check_output(['/usr/bin/objdump', '-d', '-C', '--start-address=' + str(create['address']),
                                       '--stop-address=' + str(create['address'] + create['size']), str(executable)], text=True)
    assert 'RetainedCompute::Create' in assembly
    (W / (mode + '-create.disassembly')).write_text(assembly)
    selector_names = ['kTransportPortableNormalSumShader', 'kEndpointPortableNormalSumShader',
                      'kDopriPhasePortableNormalSumShader', 'kDensePortableShader', 'kEndpointPortableShader']
    for name in selector_names:
        assert name in assembly, ('selector symbol absent from actual Create', mode, name)
    selector_records[mode] = {'symbol': create_symbols[0], 'address': symbols[create_symbols[0]]['address'],
                              'size': symbols[create_symbols[0]]['size'], 'joined_arrays': selector_names,
                              'disassembly': identity(W / (mode + '-create.disassembly'))}
for name in stock:
    assert array_records['baseline'][name]['sha256'] == hashlib.sha256(stock[name]).hexdigest()
    assert array_records['candidate'][name]['sha256'] == hashlib.sha256(candidate[name]).hexdigest()
changed = [name for name in stock if stock[name] != candidate[name]]
assert changed == [name for name in stock if name in json.loads((W / 'packet.json').read_text())['changed_arrays']]
assert len(changed) == 17
owner_array_records = {}
for mode, directory in [('baseline', B / 'src/sirius/backend/retained'), ('candidate', T / 'headers')]:
    owner_arrays = {}
    for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', (directory / 'retained_owner_kernels.h').read_text(), re.S):
        words = [int(s) for s in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count)
        owner_arrays[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    assert set(owner_arrays) == {'kTransportShader', 'kEndpointShader'}
    owner_array_records[mode], _ = elf_arrays(T / ('original_' + mode + '_controls'), owner_arrays, 'retained_owner')
contract = W / 'contract'
contract.mkdir(exist_ok=True)
reference = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse/files/out/vulkan-descriptor-reuse/audited-benefit'
for name in ['software_contracts.py', 'expected-readback-layout.json']:
    assert (contract / name).read_bytes() == (reference / name).read_bytes()
for mode, filename in [('baseline', 'baseline.json'), ('candidate', 'build-and-arrays.json')]:
    (contract / filename).write_text(json.dumps({'arrays': array_records[mode]}, indent=2) + '\n')
result = {'passed': True, 'compile_provenance': compile_records, 'link_owner': identity(W / 'link-candidate/owner.json'),
          'link_map': identity(T / 'candidate.map'), 'old_archive_compute_member_excluded': True,
          'actual_ELF_owner_arrays': owner_array_records, 'actual_ELF_arrays': array_records, 'actual_create_selectors': selector_records,
          'changed_arrays': changed, 'unchanged_array_count': 23,
          'baseline_executable': identity(T / 'original_baseline_controls'),
          'candidate_executable': identity(T / 'original_candidate_controls'),
          'original_reader': identity(contract / 'software_contracts.py'),
          'scope': 'Matched original private control producer and actual linked array/selector evidence. All seventeen affected main arrays and both unchanged optional owner arrays are bound to actual ELF bytes; execution remains for root-owned original controls. No integrated execution or benefit verdict yet.'}
(W / 'producer-acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'actual_arrays_per_ELF': 40, 'changed': changed, 'old_compute_excluded': True, 'candidate_executable': result['candidate_executable']}))
