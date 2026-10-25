"""Extend the existing ELF reader to uniquely named local-only text sections."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-accumulate-compatibility/elf_proof.py'
source = basis.read_text()
source = source.replace('if s[1] != 2:', 'if s[1] not in (2, 11):')
source = source.replace('table.append(row); self.symbols.append(row)',
    'table.append(row)\n                if s[1] == 2:\n                    self.symbols.append(row)')
before = '        functions = []\n'
assert source.count(before) == 1
source = source.replace(before, '''        # GCC emits static initialization and cold abort sections containing
        # only local functions. Use one only when its name/size is unique in
        # this exact linked ELF, then require the same section-base consistency
        # and every opcode/declared-relocation check as for public anchors.
        for symbol in obj.symbols:
            if ((symbol['info'] & 15) != 2 or not symbol['size'] or not symbol['section']
                    or bases.get(symbol['section'])):
                continue
            matches = [s for s in linked.symbols if s['name'] == symbol['name']
                       and s['section'] and s['size'] == symbol['size']]
            if len(matches) == 1:
                bases.setdefault(symbol['section'], set()).add(matches[0]['value'] - symbol['value'])
        functions = []
''')
exec(compile(source, str(basis), 'exec'))
