"""Verify this optional-ISA component without executing its arithmetic."""
from pathlib import Path
import hashlib
import json
import re
import struct
import subprocess
import sys

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BUILD = ROOT / 'bin/linux-gcc/native-cpu-v4-arithmetic'
OLD = ROOT / 'out/native-cpu-value-fact-stores'
assert WORK == ROOT / 'out/native-cpu-v4-arithmetic'
sys.path.insert(0, str(ROOT / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf, associate, pin


def doc(path):
    return json.loads(path.read_text())


def checked(path, expected):
    assert pin(path) == {key: expected[key] for key in ('bytes', 'sha256')}, path


def defined(elf, name):
    rows = [row for row in elf.symbols if row['name'] == name and row['section']]
    assert len(rows) == 1, (name, len(rows))
    return rows[0]


def compact(text):
    return re.sub(r'\s+', '', text)


derivation = doc(WORK / 'derivation.json')
compiled = doc(WORK / 'compile-result.json')
prepared = doc(OLD / 'prepared-cpp.json')
old_build = doc(OLD / 'cpp-build-result.json')
for name, expected in prepared['inputs'].items():
    checked(ROOT / name, expected)
source = (OLD / 'retained8.cpp').read_text()
checked(OLD / 'retained8.cpp', derivation['original_source'])
for name in ('optional', 'dispatch'):
    checked(WORK / (name + '.cpp'), derivation[name])
for phase in ('build', 'link'):
    owner = doc(WORK / phase / 'owner.json')
    assert owner['passed'] and owner['status'] == 'terminal'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']

signatures = re.findall(r'void (sirius_native_(?:dense|tagged)_\w+8)\((.*?)\) noexcept', source, re.S)
assert len(signatures) == 10 and [name for name, _ in signatures] == derivation['exports']
mapping = {name: name.replace('sirius_native_', 'sirius_native_v3_') for name, _ in signatures}
optional = (WORK / 'optional.cpp').read_text()
assert re.findall(r'^#define (\w+) (\w+)$', optional, re.M) == [
    (name, name.replace('sirius_native_', 'sirius_native_v4_')) for name, _ in signatures]
assert optional.splitlines()[-1] == '#include "../native-cpu-value-fact-stores/retained8.cpp"'
assert len(optional.splitlines()) == 12
assert source.count('vector_size(32)') == 3 and 'sizeof(Triple) == 160' in source
dispatch = (WORK / 'dispatch.cpp').read_text()
assert compact('bool sirius_native_v4_supported() noexcept { return __builtin_cpu_supports("x86-64-v4"); }') in compact(dispatch)
for name, parameters in signatures:
    arguments = re.findall(r'std::uint32_t\s*\*?\s+(\w+)', parameters)
    expected_arguments = ['scratch', 'tags', 'destination', 'a_index', 'b_index', 'mask'] if '_tagged_' in name else ['scratch', 'destination', 'a_index', 'b_index', 'mask']
    assert arguments == expected_arguments, (name, arguments)
    v3, v4 = mapping[name], name.replace('sirius_native_', 'sirius_native_v4_')
    for target in (v3, v4):
        assert compact('void ' + target + '(' + parameters + ') noexcept;') in compact(dispatch)
    body = re.search(r'void ' + name + r'\((.*?)\) noexcept \{(.*?)\n\}', dispatch, re.S)
    assert body and compact(body[1]) == compact(parameters)
    args = ','.join(arguments)
    assert compact(body[2]) == compact('if(sirius_native_v4_supported()) ' + v4 + '(' + args + '); else ' + v3 + '(' + args + ');')
assert dispatch.count('extern "C"') == 31

old_flags = list(prepared['compilations'][0]['argv'])
old_flags = old_flags[:old_flags.index('-c')]
assert len(compiled['commands']) == 2
for label, row in zip(('optional', 'dispatch'), compiled['commands']):
    flags = list(old_flags)
    flags[flags.index('-march=x86-64-v3')] = '-march=x86-64-v4' if label == 'optional' else '-march=x86-64'
    if label == 'optional':
        flags.append('-mprefer-vector-width=256')
    assert row['exit'] == 0 and row['argv'] == flags + ['-c', str(WORK / (label + '.cpp')), '-o', str(BUILD / (label + '.o'))]
    assert not any(flag.startswith('-flto') for flag in flags)
for label, expected in compiled['objects'].items():
    checked(BUILD / (label + '.o'), expected)
original = ROOT / compiled['original_v3']['path']
checked(original, compiled['original_v3'])
a, b = Elf(original), Elf(BUILD / 'v3.o')
inverse = {new: old for old, new in mapping.items()}
assert sorted(a.descriptor(row, {}) for row in a.symbols) == sorted(b.descriptor(row, inverse) for row in b.symbols)
allocated = {}
for index, section in enumerate(a.sections):
    if not section[2] & 2:
        continue
    candidates = [i for i, name in enumerate(b.names) if name == a.names[index]]
    assert len(candidates) == 1
    other = candidates[0]
    assert tuple(section[i] for i in (1, 2, 3, 5, 8, 9)) == tuple(b.sections[other][i] for i in (1, 2, 3, 5, 8, 9))
    assert a.bytes(index) == b.bytes(other)
    allocated[a.names[index]] = dict(bytes=section[5], sha256=hashlib.sha256(a.bytes(index)).hexdigest())
assert set(allocated) == {b.names[i] for i, section in enumerate(b.sections) if section[2] & 2}
relocs = lambda elf, inv: sorted((elf.section_name(section), offset, kind, addend, elf.descriptor(target, inv)) for section, offset, kind, addend, target in elf.relocations)
assert relocs(a, {}) == relocs(b, inverse)


def isa(path, optional_isa=False, baseline=False):
    elf = Elf(path)
    assert not any(section[5] and (name.startswith(('.init', '.preinit', '.ctors')))
                   for name, section in zip(elf.names, elf.sections))
    assert not any((row['info'] & 15) == 10 or row['name'].startswith(('_GLOBAL__sub_I_', '__static_initialization')) for row in elf.symbols)
    raw = subprocess.check_output(['objdump', '-d', '--insn-width=16', str(path)], text=True)
    mnemonics, evex = set(), 0
    for line in raw.splitlines():
        fields = line.split('\t')
        if len(fields) < 3 or not re.match(r'^\s*[0-9a-f]+:$', fields[0]):
            continue
        encoded, text = fields[1].split(), fields[2].strip()
        assert encoded and text
        evex += encoded[0] == '62'
        mnemonics.add(text.split()[0])
    forbidden = {m for m in mnemonics if m.startswith(('vfm', 'vfnm', 'f', 'aes', 'vaes', 'pclmul', 'vpclmul', 'sha', 'vgf2p8', 'vpdp', 'vpmadd52', 'tile', 'tdp')) or m in ('ldmxcsr', 'stmxcsr', 'vldmxcsr', 'vstmxcsr', 'rdrand', 'rdseed', 'adcx', 'adox', 'rdpid', 'rdtscp', 'prefetchw', 'prefetchwt1', 'clflushopt', 'clwb', 'clzero', 'movdiri', 'movdir64b', 'serialize', 'pconfig')}
    assert not forbidden, (path, forbidden)
    assert not re.search(r'%(?:mm[0-7]|st|tmm\d+)\b', raw)
    upper = sorted(set(re.findall(r'%[xy]mm(?:1[6-9]|2\d|3[01])\b', raw)))
    if optional_isa:
        assert evex and upper
    else:
        assert not evex and not upper and not re.search(r'%(?:zmm\d+|k[0-7])\b', raw)
    if baseline:
        assert mnemonics <= {'endbr64', 'mov', 'movzbl', 'shr', 'and', 'ret', 'cs', 'testb', 'je', 'jmp', 'data16'}
    return dict(object=pin(path), EVEX_instructions=evex, upper_vector_registers=upper,
                mnemonics=sorted(mnemonics), no_initializer_or_IFUNC=True)


partitions = {label: isa(BUILD / (label + '.o'), label == 'optional', label == 'dispatch')
              for label in ('v3', 'optional', 'dispatch')}
items = [dict(path=str(BUILD / (label + '.o')), **compiled['objects'][label]) for label in ('v3', 'optional')]
joins = associate(BUILD / 'abi-check', items)
old_object = dict(path=str(original), **pin(original))
original_join = associate(ROOT / old_build['executable']['path'], [old_object])
linked, dispatcher = Elf(BUILD / 'abi-check'), Elf(BUILD / 'dispatch.o')
cpu = defined(linked, '__cpu_features2')
dispatch_joins = []
for symbol in dispatcher.symbols:
    if (symbol['info'] & 15) != 2 or not symbol['size'] or not symbol['section']:
        continue
    name = symbol['name']
    final_symbol = defined(linked, name)
    assert final_symbol['size'] == symbol['size']
    before, after = bytearray(dispatcher.function_bytes(symbol)), bytearray(linked.function_bytes(final_symbol))
    assert before[:7] == b'\xf3\x0f\x1e\xfa\x48\x8b\x05'
    assert after[:7] == b'\xf3\x0f\x1e\xfa\x48\x8d\x05'
    assert final_symbol['value'] + 11 + struct.unpack_from('<i', after, 7)[0] == cpu['value']
    after[5] = before[5]  # The one checked GOTPCRELX mov-to-lea relaxation.
    relocations = [(offset - symbol['value'], kind, addend, target['name'])
                   for section, offset, kind, addend, target in dispatcher.relocations
                   if section == symbol['section'] and 0 <= offset - symbol['value'] < symbol['size']]
    expected = [(7, 42, -4, '__cpu_features2')]
    targets = []
    if name == 'sirius_native_v4_supported':
        assert before[11:] == bytes.fromhex('0fb64008c0e80283e001c3')
    else:
        assert name in mapping and len(before) == 37
        assert before[11:18] == bytes.fromhex('f6400804740fe9')
        assert before[22:32] == bytes.fromhex('662e0f1f840000000000') and before[32] == 0xe9
        for offset, target in ((18, name.replace('sirius_native_', 'sirius_native_v4_')), (33, mapping[name])):
            expected.append((offset, 4, -4, target))
            value = final_symbol['value'] + offset + 4 + struct.unpack_from('<i', after, offset)[0]
            assert value == defined(linked, target)['value'], (name, target)
            targets.append(dict(relocation_offset=offset, target=target, linked_address=value))
    assert sorted(relocations) == sorted(expected)
    for offset, _kind, _addend, _target in relocations:
        before[offset:offset + 4] = after[offset:offset + 4] = b'\0' * 4
    assert before == after, name
    dispatch_joins.append(dict(function=name, linked_address=final_symbol['value'],
        bytes=symbol['size'], cpu_features2_address=cpu['value'],
        checked_GOTPCRELX_mov_to_lea=True, false_branch_v3_true_branch_v4=name in mapping,
        targets=targets))
assert len(dispatch_joins) == 11
notes = subprocess.check_output(['readelf', '-n', str(BUILD / 'abi-check')], text=True)
assert re.findall(r'x86 ISA needed:\s*([^\n]+)', notes) == ['x86-64-baseline']
assert not re.search(r'x86-64-v[234]', notes)
result = dict(passed=True, source_exports_verified=10, unchanged_original_v3=pin(original),
    v3_inverse_symbols_relocations_and_allocated_bytes_equal=True, allocated_sections=allocated,
    strict_compile_commands_verified=2, ISA_partitions=partitions,
    original_function_joins=original_join, new_math_function_joins=joins,
    dispatcher_function_joins=dispatch_joins, linked_ELF_notes=notes,
    verifier=pin(Path(__file__)),
    limits='Source, allocated bytes, symbols, declared relocations, partitioned ISA and exact observed dispatcher control-flow only. '
           'GOTPCRELX mov-to-lea is checked against the actual __cpu_features2 symbol; no arbitrary opcode relaxation is masked. '
           'GCC builtin CPU initialization/OS detection is an external compiler-runtime contract. '
           'This does not establish executed v4 PCs, execution on older hardware, arithmetic correctness, runtime gain, provider adoption, frame completion or scientific/release qualification.')
with (WORK / 'code-proof.json').open('x') as stream:
    stream.write(json.dumps(result, indent=2) + '\n')
print(json.dumps(dict(passed=True, source_exports=10, dispatcher_joins=len(dispatch_joins),
    original_function_joins=len(original_join['objects'][0]['functions']),
    new_math_function_joins=sum(len(row['functions']) for row in joins['objects']))))
