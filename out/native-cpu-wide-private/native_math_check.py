"""Reuse the retained full-FDE reader for each new regcall value adapter."""
from pathlib import Path
import hashlib
import re

R = Path.cwd()
BASE = R / 'out/native-cpu-canonical-closed-remaining/native_check.py'
assert hashlib.sha256(BASE.read_bytes()).hexdigest() == '99852530e00968fc81f577774e224791fcba4c197e8b8a436ac41914f0c5303e'
source = BASE.read_text()
before = "    assert len(functions) == 1\n"
assert source.count(before) == 1
source = source.replace(before, "    functions = [row for row in functions if row['name'] == symbol]\n" + before, 1)
before = "    instructions = [(int(row[0], 16) - start, row[1])"
assert source.count(before) == 1
source = source.replace(before,
    "    pattern = r'^[0-9a-f]+ <' + re.escape(symbol) + r'>:\\n(.*?)(?=^[0-9a-f]+ <|\\Z)'\n"
    "    section = re.search(pattern, assembly, re.M | re.S)\n"
    "    assert section is not None\n"
    "    assembly = section[1]\n" + before, 1)
source = source.replace('sole_hidden_C_ABI_entry=symbol', 'hidden_regcall_value_entry=symbol')
source = source.replace('sole native entry', 'selected value entry')
namespace = {'__file__': str(BASE), '__name__': 'retained_native_FDE_reader'}
exec(compile(source, str(BASE), 'exec'), namespace)
check = namespace['check']


def all_adapters(path, names):
    elf = namespace['Elf'](path)
    functions = [row for row in elf.symbols if row['section'] and row['info'] & 15 == 2 and row['size']]
    assert len(functions) == len(names) == 9
    assert {row['name'] for row in functions} == set(names)
    return {name: check(path, name) for name in names}
