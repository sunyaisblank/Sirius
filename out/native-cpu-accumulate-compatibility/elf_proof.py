"""Read-only ELF64/x86_64 rename and relocation-masked function association."""
from pathlib import Path
import hashlib
import struct


def pin(path):
    with Path(path).open('rb') as stream:
        return {'bytes': Path(path).stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


class Elf:
    def __init__(self, path):
        self.data = Path(path).read_bytes()
        h = struct.unpack_from('<16sHHIQQQIHHHHHH', self.data)
        assert h[0][:7] == b'\x7fELF\x02\x01\x01' and h[2] == 62 and h[11] == 64
        self.kind = h[1]
        self.sections = [struct.unpack_from('<IIQQQQIIQQ', self.data, h[6] + i * 64) for i in range(h[12])]
        names = self.bytes(h[13])
        self.names = [self.string(names, s[0]) for s in self.sections]
        self.tables, self.symbols = {}, []
        for i, s in enumerate(self.sections):
            if s[1] != 2:
                continue
            assert s[9] == 24
            strings = self.bytes(s[6]); table = []
            for j in range(0, s[5], 24):
                name, info, other, section, value, size = struct.unpack_from('<IBBHQQ', self.bytes(i), j)
                row = dict(name=self.string(strings, name), info=info, other=other,
                           section=section, value=value, size=size)
                table.append(row); self.symbols.append(row)
            self.tables[i] = table
        self.relocations = []
        for i, s in enumerate(self.sections):
            assert s[1] != 9, 'REL without explicit addends is unsupported'
            if s[1] == 4:
                assert s[9] == 24
                for j in range(0, s[5], 24):
                    offset, info, addend = struct.unpack_from('<QQq', self.bytes(i), j)
                    self.relocations.append((s[7], offset, info & 0xffffffff, addend,
                                             self.tables[s[6]][info >> 32]))

    @staticmethod
    def string(data, offset):
        return data[offset:data.index(b'\0', offset)].decode()

    def bytes(self, index):
        s = self.sections[index]
        return b'' if s[1] == 8 else self.data[s[4]:s[4] + s[5]]

    def section_name(self, index):
        return self.names[index] if index < len(self.names) else 'SPECIAL:' + str(index)

    def function_bytes(self, row):
        s = self.sections[row['section']]
        first = s[4] + row['value'] - s[3]
        return self.data[first:first + row['size']]

    def descriptor(self, row, inverse):
        return (inverse.get(row['name'], row['name']), row['info'], row['other'],
                self.section_name(row['section']), row['value'], row['size'])


def rename_proof(original, renamed, mapping):
    a, b = Elf(original), Elf(renamed); assert a.kind == b.kind == 1
    inverse = {v: k for k, v in mapping.items()}; assert len(inverse) == len(mapping) == 2
    def definitions(elf):
        return {s['name'] for s in elf.symbols if s['section'] and s['info'] >> 4 in (1, 2)}
    assert definitions(a) == set(mapping) and definitions(b) == set(inverse)
    def undefined(elf):
        return {s['name'] for s in elf.symbols if not s['section'] and s['name']}
    assert undefined(a) == undefined(b) == {'memset', 'sirius_simd_fmaxf'}
    assert sorted(a.descriptor(s, {}) for s in a.symbols) == sorted(b.descriptor(s, inverse) for s in b.symbols)
    allocated = {}
    for i, s in enumerate(a.sections):
        if not s[2] & 2:
            continue
        matches = [j for j, t in enumerate(b.sections) if b.names[j] == a.names[i]]
        assert len(matches) == 1
        j = matches[0]; t = b.sections[j]
        assert (s[1], s[2], s[3], s[5], s[8], s[9]) == (t[1], t[2], t[3], t[5], t[8], t[9])
        assert a.bytes(i) == b.bytes(j)
        allocated[a.names[i]] = {'bytes': s[5], 'type': s[1], 'sha256': hashlib.sha256(a.bytes(i)).hexdigest()}
    assert {a.names[i] for i, s in enumerate(a.sections) if s[2] & 2} == {
        b.names[i] for i, s in enumerate(b.sections) if s[2] & 2}
    def relocs(elf, inverse):
        return sorted((elf.section_name(section), offset, kind, addend, elf.descriptor(target, inverse))
                      for section, offset, kind, addend, target in elf.relocations)
    assert relocs(a, {}) == relocs(b, inverse)
    return {'passed': True, 'original': pin(original), 'renamed': pin(renamed), 'mapping': mapping,
            'all_allocated_sections_identical': allocated, 'all_symbols_inverse_normalized_equal': True,
            'all_relocations_inverse_normalized_equal': True, 'undefined_symbols_unchanged': sorted(undefined(a))}


def associate(executable, objects):
    linked = Elf(executable); rows = []
    for item in objects:
        path = Path(item['path']); assert pin(path) == {k: item[k] for k in ('bytes', 'sha256')}
        obj = Elf(path); assert obj.kind == 1
        bases = {}
        for symbol in obj.symbols:
            if (symbol['info'] & 15) != 2 or not symbol['size'] or not symbol['section']:
                continue
            if symbol['info'] >> 4 not in (1, 2):
                continue
            matches = [s for s in linked.symbols if s['name'] == symbol['name'] and s['section']]
            assert len(matches) == 1 and matches[0]['size'] == symbol['size']
            bases.setdefault(symbol['section'], set()).add(matches[0]['value'] - symbol['value'])
        functions = []
        for symbol in obj.symbols:
            if (symbol['info'] & 15) != 2 or not symbol['size'] or not symbol['section']:
                continue
            section_bases = bases.get(symbol['section'], set()); assert len(section_bases) == 1
            address = next(iter(section_bases)) + symbol['value']
            matches = [s for s in linked.symbols if s['name'] == symbol['name'] and
                       s['section'] and s['value'] == address and s['size'] == symbol['size']]
            assert len(matches) == 1, (path, symbol['name'], address)
            original, final = bytearray(obj.function_bytes(symbol)), bytearray(linked.function_bytes(matches[0]))
            masked = []
            for section, location, kind, _addend, _target in obj.relocations:
                offset = location - symbol['value']
                if section == symbol['section'] and 0 <= offset < symbol['size']:
                    assert kind in (2, 4, 42) and offset + 4 <= symbol['size']
                    original[offset:offset + 4] = final[offset:offset + 4] = b'\0' * 4
                    masked.append({'offset': offset, 'ELF_type': kind, 'bytes': 4})
            assert original == final, (path, symbol['name'])
            functions.append({'symbol': symbol['name'], 'bytes': symbol['size'], 'linked_address': address,
                              'masked_relocations': masked})
        assert functions
        rows.append({'object': item, 'all_retained_function_opcodes_match': True, 'functions': functions})
    return {'passed': True, 'executable': {'path': str(executable), **pin(executable)}, 'objects': rows,
            'scope': 'Every retained function opcode in exact renamed OLD/new R Transport objects, masking only declared4B PC32/PLT32/REX_GOTPCRELX operands. Linker relaxation disabled. This alone is not relocation-target semantics or numerical acceptance.'}
