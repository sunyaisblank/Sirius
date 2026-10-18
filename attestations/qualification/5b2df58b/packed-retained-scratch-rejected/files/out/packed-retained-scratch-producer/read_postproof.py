"""Read existing matched producer data only; no imports of project code or tools."""
from pathlib import Path
import hashlib
import json
import re
import shlex
import struct

W = Path(__file__).resolve().parent
R = W.parents[1]
S = R / 'out/packed-retained-scratch-controls'
B = R / 'bin/linux-gcc'
T = B / 'packed-retained-scratch-controls'
digest_cache = {}

def document(path):
    return json.loads(path.read_text())

def identity(path):
    resolved = str(path.resolve())
    if resolved not in digest_cache:
        data = path.read_bytes()
        digest_cache[resolved] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    return digest_cache[resolved]

def pin(path):
    return {'path': str(path.relative_to(R)), **identity(path)}

def same(path, expected):
    assert identity(path) == {k: expected[k] for k in ('bytes', 'sha256')}, str(path)
    if 'resolved_path' in expected:
        assert str(path.resolve()) == expected['resolved_path'], str(path)

def arrays(path):
    result = {}
    pattern = r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};'
    for count, name, body in re.findall(pattern, path.read_text(), re.S):
        words = [int(word) for word in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count) and name not in result
        result[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    return result

def elf_symbols(path):
    data = path.read_bytes()
    header = struct.unpack_from('<16sHHIQQQIHHHHHH', data)
    assert header[0][:6] == b'\x7fELF\x02\x01' and header[2] == 62 and header[11] == 64
    assert header[6] + header[11] * header[12] <= len(data)
    sections = [struct.unpack_from('<IIQQQQIIQQ', data, header[6] + i * header[11])
                for i in range(header[12])]
    result = {}
    for section in sections:
        if section[1] != 2:
            continue
        assert section[9] == 24 and section[4] + section[5] <= len(data)
        strings = sections[section[6]]
        names = data[strings[4]:strings[4] + strings[5]]
        for off in range(section[4], section[4] + section[5], 24):
            name, info, other, index, value, size = struct.unpack_from('<IBBHQQ', data, off)
            assert name < len(names)
            end = names.find(b'\0', name)
            assert end >= name
            label = names[name:end].decode()
            if not label or not 0 < index < len(sections):
                continue
            sec = sections[index]
            raw = None
            if sec[1] != 8 and value >= sec[3] and value + size <= sec[3] + sec[5]:
                begin = sec[4] + value - sec[3]
                assert begin + size <= len(data)
                raw = data[begin:begin + size]
            record = {'address': value, 'size': size, 'binding': info >> 4,
                      'type': info & 15, 'section_flags': sec[2], 'data': raw}
            if label.startswith(('_ZN6sirius7backend16retained_program',
                                 '_ZN6sirius7backend14retained_owner',
                                 '_ZN6sirius7backend15RetainedCompute6Create')):
                assert label not in result, ('duplicate required symbol', label)
            result[label] = record
    return result

def check_arrays(symbols, namespace, expected):
    records = {}
    prefix = '_ZN6sirius7backend' + str(len(namespace)) + namespace
    for name, data in expected.items():
        symbol = symbols[prefix + str(len(name)) + name + 'E']
        assert symbol['type'] == 1 and symbol['section_flags'] & 2
        assert symbol['size'] == len(data) and symbol['data'] == data, name
        records[name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
                         'address': symbol['address'], 'binding': symbol['binding']}
    return records

preflight = document(W / 'root-preflight.json')
action = document(W / 'actions.json')['producer']
owner = document(W / 'producer/owner.json')
before = document(W / 'producer/inputs-before.json')
after = document(W / 'producer/inputs-after.json')
assert before == after and len(before) == 5522 and len(preflight['files']) == 5521
assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert owner['argv'] == action['argv'] == preflight['argv']
assert owner['limits'] == action['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
assert owner['stop_reason'] is None and not owner['cleanup_errors']
assert owner['remaining_owned_processes'] == [] and owner['provider_mappings'] == []
assert owner['source'] == {'revision': preflight['source_revision'], 'status': ''}
normalized = {str((R / name).resolve()): row for name, row in before.items()}
actual_paths = {str(R / name) for name in before}
assert actual_paths == set(preflight['files']) | {str(W / 'root-preflight.json')}
for name, row in before.items():
    same(R / name, row)
for name, row in preflight['files'].items():
    same(Path(name), row)

# Fresh checks are scoped strictly to recorded births, with PID reuse distinguished.
assert Path('/proc/sys/kernel/random/boot_id').read_text().strip() == owner['boot_id']
births = {(int(row['pid']), str(row['saved_start_ticks'])) for row in owner['observed_births']}
assert all(row['absent'] for row in owner['observed_births'])
births.add((int(owner['controller']['pid']), str(owner['controller']['start_ticks'])))
for pid, saved in births:
    try:
        stat = Path('/proc') / str(pid) / 'stat'
        fields = stat.read_text().rpartition(') ')[2].split()
        assert fields[19] != saved, ('recorded birth remains live', pid, saved)
    except FileNotFoundError:
        pass

samples = [json.loads(line) for line in (W / 'producer/samples.jsonl').read_text().splitlines()]
assert samples
sample_keys = sorted(samples[0])
assert max(row['rss_kib'] for row in samples) == owner['peak_sampled_rss_kib']
assert min(row['available_kib'] for row in samples) == owner['minimum_sampled_available_kib']
assert owner['peak_sampled_rss_kib'] <= owner['limits']['rss_kib']
assert owner['minimum_sampled_available_kib'] >= owner['limits']['reserve_kib']
assert owner['wall_seconds'] < owner['limits']['seconds']
recipes = document(S / 'recipes.pending.json')
commands = document(T / 'producer-command-results.json')
labels = recipes['producer_operation_order']
assert commands['passed'] and len(commands['commands']) == len(labels) == 6
for row, label in zip(commands['commands'], labels):
    assert row == {'label': label, 'argv': recipes['producer_actions'][label]['argv'],
                   'started': True, 'exit': 0}
printed = [json.loads(line) for line in (W / 'producer/stdout.log').read_text().splitlines()
           if line.startswith('{"starting":')]
assert printed == [{'starting': label, 'argv': recipes['producer_actions'][label]['argv']}
                   for label in labels]
assert len(commands['products']) == 12 and len(commands['resources']) == 23
for name, row in commands['products'].items():
    same(R / name, row)
resources = {str(path.relative_to(R)): identity(path)
             for path in sorted((T / 'resources').rglob('*')) if path.is_file()}
assert resources == commands['resources']
for name, row in resources.items():
    original = B / 'tests/backend/resources' / Path(name).relative_to(T.relative_to(R) / 'resources')
    same(original, row)
    assert str(original.resolve()) in normalized

dependencies = {}
for arm in ('baseline', 'candidate'):
    for kind in ('compute', 'observer'):
        depfile = T / (arm + '-' + kind + '.d')
        target, body = depfile.read_text().replace('\\\n', '').split(':', 1)
        assert target.strip() == str(T / (arm + '-' + kind + '.o'))
        paths = {str((R / value).resolve()) for value in shlex.split(body)}
        header = str((S / 'headers' / arm / 'retained_kernels.h').resolve())
        other = str((S / 'headers' / ('candidate' if arm == 'baseline' else 'baseline') / 'retained_kernels.h').resolve())
        assert header in paths and other not in paths
        for name in paths:
            assert name in normalized, ('unsealed actual compiler dependency', arm, kind, name)
            same(Path(name), normalized[name])
        dependencies[arm + '-' + kind] = {'depfile': pin(depfile), 'count': len(paths),
                                         'correct_private_header': pin(Path(header)),
                                         'all_dependencies_sealed_and_current': True}

stock = arrays(S / 'headers/baseline/retained_kernels.h')
candidate = arrays(S / 'headers/candidate/retained_kernels.h')
assert len(stock) == len(candidate) == 40
changed = sorted(name for name in stock if stock[name] != candidate[name])
assert len(changed) == 25
programs = sorted(name for name in stock if name.endswith('Program'))
assert len(programs) == 7 and all(stock[name] == candidate[name] for name in programs)
owner_expected = arrays(B / 'src/sirius/backend/retained/retained_owner_kernels.h')
assert len(owner_expected) == 2
constant_pattern = r'inline constexpr (std::size_t|bool) (\w+) = ([^;]+);'
constants = re.findall(constant_pattern, (S / 'headers/baseline/retained_kernels.h').read_text())
assert constants and constants == re.findall(constant_pattern, (S / 'headers/candidate/retained_kernels.h').read_text())

device_path = B / 'src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o'
assert str(device_path.resolve()) in normalized
device_symbols = elf_symbols(device_path)
strong_device = sorted(name for name, row in device_symbols.items()
                       if row['binding'] == 1 and row['type'] == 2)
assert strong_device
elf_arrays = {}
elf_owner = {}
creates = {}
map_records = {}
linked_device_counts = {}
for arm, expected in [('baseline', stock), ('candidate', candidate)]:
    executable = T / ('original_' + arm + '_controls')
    symbols = elf_symbols(executable)
    elf_arrays[arm] = check_arrays(symbols, 'retained_program', expected)
    elf_owner[arm] = check_arrays(symbols, 'retained_owner', owner_expected)
    for name in strong_device:
        assert symbols[name]['binding'] == 1 and symbols[name]['type'] == 2 and symbols[name]['size'] > 0
    linked_device_counts[arm] = len([name for name in strong_device
                                   if symbols[name]['binding'] == 1 and symbols[name]['type'] == 2])
    assert linked_device_counts[arm] == len(strong_device)
    create_names = [name for name in symbols if name.startswith('_ZN6sirius7backend15RetainedCompute6Create') and '.' not in name]
    assert len(create_names) == 1
    create = symbols[create_names[0]]
    assert create['type'] == 2 and create['binding'] == 1 and create['size'] > 0
    creates[arm] = {'symbol': create_names[0], 'address': create['address'], 'size': create['size'],
                    'binding': create['binding'], 'executable': pin(executable)}
    linkmap = (T / (arm + '.map')).read_text()
    for member in ('retained_compute.cpp.o', 'vulkan_device.cpp.o'):
        assert not re.search(r'libsirius_backend\.a\([^)]*' + re.escape(member) + r'\)', linkmap)
    for name in (T / (arm + '-compute.o'), T / (arm + '-observer.o'), device_path):
        assert str(name) in linkmap
    assert '--whole-archive' not in recipes['producer_actions']['link-' + arm]['argv']
    map_records[arm] = pin(T / (arm + '.map'))

posttools = []
for arm, create in creates.items():
    executable = str(T / ('original_' + arm + '_controls'))
    posttools += [
        {'label': arm + '-create', 'argv': ['/usr/bin/objdump', '-d', '-C',
         '--start-address=' + str(create['address']), '--stop-address=' + str(create['address'] + create['size']), executable]},
        {'label': arm + '-symbols', 'argv': ['/usr/bin/nm', '-C', '--defined-only', executable]},
        {'label': arm + '-readelf', 'argv': ['/usr/bin/readelf', '-Ws', executable]}]
(W / 'posttools.pending.json').write_text(json.dumps({'schema': 'sirius-packed-scratch-posttools-pending-v1',
    'executed': False, 'root_pin_and_owned_batch_required': True, 'commands': posttools,
    'required_check': 'Actual Create selector references to respective ELF array addresses; tool output not yet available.'}, indent=2) + '\n')

result = {'schema': 'sirius-packed-scratch-current-producer-data-postproof-v1', 'passed': True,
    'source_revision': preflight['source_revision'], 'method': 'Pure stdlib data reader; no project imports, tools, executables or provider/test launches.',
    'reader': pin(Path(__file__)), 'owner': pin(W / 'producer/owner.json'),
    'input_maps': {'before': pin(W / 'producer/inputs-before.json'), 'after': pin(W / 'producer/inputs-after.json'),
                   'rows_each': len(before), 'fixed_rich_rows': len(preflight['files']), 'all_current': True},
    'ownership': {'same_boot_recorded_births_fresh_retired': len(births), 'observed_births': len(owner['observed_births']),
                  'wall_seconds': owner['wall_seconds'], 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'],
                  'minimum_sampled_available_kib': owner['minimum_sampled_available_kib'],
                  'sample_count': len(samples), 'sample_keys': sample_keys, 'stop_reason': None, 'cleanup_errors': [], 'remaining_owned_processes': []},
    'commands': {'receipt': pin(T / 'producer-command-results.json'), 'all_six_argv_and_exits_and_stdout_match': True,
                 'stdout': pin(W / 'producer/stdout.log'), 'stderr': pin(W / 'producer/stderr.log'),
                 'samples': pin(W / 'producer/samples.jsonl'),
                 'products': commands['products'], 'resources_exact_current_sealed_copies': len(resources)},
    'compile_dependencies': dependencies, 'link_maps': map_records,
    'default_archive_compute_and_device_members_excluded': True,
    'actual_ELF_arrays': elf_arrays, 'actual_ELF_owner_arrays': elf_owner,
    'changed_arrays': changed, 'unchanged_main_array_count': 15, 'seven_tables': programs,
    'public_constants': {'header_declarations_byte_exact': constants,
                         'scope': 'Source constants are constexpr and may be folded; no claim to independently recover every folded constant from machine code.'},
    'strong_device': {'object': pin(device_path), 'actual_count': len(strong_device), 'symbols': strong_device,
                       'actual_linked_counts': linked_device_counts,
                       'same_defined_global_function_binding_in_both_ELFs': True},
    'create_functions': creates, 'posttools_pending': pin(W / 'posttools.pending.json'),
    'create_selector_disassembly_proved': False, 'runtime_admitted': False,
    'scope': 'Current matched C++ producer/input/ELF data proof only. Root review and pinned GNU selector batch still required before runtime. No original control, ABBA, numerical, driver-cost, adoption or qualification result.'}
(W / 'producer-postproof.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'input_rows': len(before), 'six_commands': 6,
                  'arrays_each': 40, 'owner_arrays_each': 2, 'changed': len(changed),
                  'strong_device_functions_actual': len(strong_device),
                  'recorded_births_retired': len(births)}))
