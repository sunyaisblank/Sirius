"""Bind linked native function bytes to sealed objects, allowing ELF relocations."""
from pathlib import Path
import hashlib
import json
import re
import struct
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


class Elf:
    def __init__(self, path):
        self.data = path.read_bytes()
        header = struct.unpack_from('<16sHHIQQQIHHHHHH', self.data)
        assert header[0][:7] == b'\x7fELF\x02\x01\x01' and header[2] == 62
        assert header[11] == 64
        self.sections = [struct.unpack_from('<IIQQQQIIQQ', self.data, header[6] + i * 64)
                         for i in range(header[12])]
        names = self.section_data(header[13])
        self.names = [self.string(names, section[0]) for section in self.sections]
        self.symbols = []
        for index, section in enumerate(self.sections):
            if section[1] != 2:
                continue
            strings = self.section_data(section[6])
            for offset in range(0, section[5], 24):
                name, info, other, shndx, value, size = struct.unpack_from(
                    '<IBBHQQ', self.section_data(index), offset)
                self.symbols.append({'name': self.string(strings, name), 'type': info & 15,
                                     'section': shndx, 'value': value, 'size': size})

    @staticmethod
    def string(data, offset):
        return data[offset:data.index(b'\0', offset)].decode()

    def section_data(self, index):
        section = self.sections[index]
        return self.data[section[4]:section[4] + section[5]]

    def function_data(self, symbol):
        section = self.sections[symbol['section']]
        begin = section[4] + symbol['value'] - section[3]
        return self.data[begin:begin + symbol['size']]

    def named(self, name):
        matches = [symbol for symbol in self.symbols if symbol['name'] == name and symbol['section']]
        assert len(matches) == 1, (name, len(matches))
        return matches[0]


build = json.loads((P / 'provider-build-result.json').read_text())
library = R / build['library']
assert pin(library) == build['library_identity']
linked = Elf(library)
rows = []
for item in build['native_objects']:
    path = R / item['path']
    assert pin(path) == {key: item[key] for key in ('bytes', 'sha256')}
    obj = Elf(path)
    relocations = {}
    for index, section in enumerate(obj.sections):
        if section[1] != 4:
            continue
        for offset in range(0, section[5], 24):
            location, info, _addend = struct.unpack_from('<QQq', obj.section_data(index), offset)
            relocations.setdefault(section[7], []).append((location, info & 0xffffffff))
    functions = []
    omitted_noop_initializers = []
    for symbol in obj.symbols:
        if symbol['type'] != 2 or not symbol['size'] or not symbol['section']:
            continue
        if symbol['name'].startswith('_GLOBAL__sub_I_'):
            assert obj.function_data(symbol) == b'\xc3'
            assert not [entry for entry in linked.symbols if entry['name'] == symbol['name']]
            omitted_noop_initializers.append(symbol['name'])
            continue
        final = linked.named(symbol['name'])
        assert final['size'] == symbol['size'], symbol['name']
        original_bytes = bytearray(obj.function_data(symbol))
        final_bytes = bytearray(linked.function_data(final))
        allowed_relocations = []
        for location, kind in relocations.get(symbol['section'], []):
            offset = location - symbol['value']
            if 0 <= offset < symbol['size']:
                assert kind in {2, 4, 42} and offset + 4 <= symbol['size'], (symbol['name'], kind)
                original_bytes[offset:offset+4] = b'\0' * 4
                final_bytes[offset:offset+4] = b'\0' * 4
                allowed_relocations.append({'offset': offset, 'ELF_type': kind, 'bytes': 4})
        assert original_bytes == final_bytes, symbol['name']
        functions.append({'symbol': symbol['name'], 'bytes': symbol['size'],
                          'linked_address': final['value'], 'relocations': allowed_relocations})
    assert functions
    rows.append({'object': item, 'all_retained_function_opcodes_match': True,
                 'omitted_private_one_byte_RET_initializers': omitted_noop_initializers,
                 'functions': functions})

disassembly = subprocess.check_output(['objdump', '-d', '--disassemble=cs_exec_fn', str(library)], text=True)
for slug in ('transport', 'endpoint', 'dopri_phase'):
    assert '<sirius_native_' + slug + '_try_packet>' in disassembly
assert '<util_fpstate_get>' in disassembly and disassembly.count('<util_fpstate_set>') == 2
assert '$0xffff003f' in disassembly and '$0x1f80' in disassembly
for name, instruction in [('util_fpstate_get', 'stmxcsr'), ('util_fpstate_set', 'ldmxcsr')]:
    helper = subprocess.check_output(['objdump', '-d', '--disassemble=' + name, str(library)], text=True)
    assert instruction in helper
    disassembly += helper
(P / 'linked-callback.disassembly.txt').write_text(disassembly)
all_native_disassembly = subprocess.check_output(['objdump', '-d', str(library),
    '--start-address=' + str(min(function['linked_address'] for row in rows for function in row['functions']
                                if '_GLOBAL__sub_I_' not in function['symbol'])),
    '--stop-address=' + str(max(function['linked_address'] + function['bytes']
                               for row in rows for function in row['functions']
                               if '_GLOBAL__sub_I_' not in function['symbol']))], text=True)
assert not re.search(r'\b(?:v?fm(?:add|sub)|v?fnm(?:add|sub))\w*\b', all_native_disassembly)
(P / 'linked-native.disassembly.txt').write_text(all_native_disassembly)
resolutions = subprocess.check_output(['ldd', '-r', str(library)], text=True, stderr=subprocess.STDOUT)
assert 'undefined symbol' not in resolutions and 'not found' not in resolutions
(P / 'linked-runtime-resolution.txt').write_text(resolutions)
batch = json.loads((W / 'batch-generation.json').read_text())
for stage in ('Transport', 'Endpoint', 'DopriPhase'):
    module = (R / batch['modules'][stage]['SPIRV']['module']).read_bytes()
    assert linked.data.count(module) == 1, stage
(P / 'linked-native-proof.json').write_text(json.dumps({
    'passed': True, 'library': build['library'], 'library_identity': pin(library),
    'objects': rows, 'callback_direct_calls_all3': True,
    'callback_saved_and_restored_MXCSR_instructions_present': True,
    'native_FMA_mnemonics': 0, 'all3_complete_SPIRV_spans_present_once': True,
    'linked_runtime_symbols_resolve': True,
    'scope': 'Static ELF/code association only: all retained sealed object function opcodes equal linked code after masking only declared4B PC32/PLT32/REX_GOTPCRELX relocation operands. Three omitted private one-byteRET initializers have no work. Does not alone validate relocation semantics, dynamic native admission, libm/FP equivalence, numerical or benefit/frame acceptance. Callback includes direct native calls and ordinary fallback under existing IEEE32 guard.',
}, indent=2) + '\n')
print(json.dumps({'linked_functions': sum(len(row['functions']) for row in rows),
                  'function_opcodes_match': True, 'runtime_symbols_resolve': True}))
