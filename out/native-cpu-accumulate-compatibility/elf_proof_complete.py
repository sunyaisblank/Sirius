"""Read dynamic symbol tables without treating them as duplicate static symbols."""
from pathlib import Path
import hashlib

basis = Path(__file__).resolve().parent / 'elf_proof.py'
source = basis.read_text()
before = "if s[1] != 2:"
assert source.count(before) == 1
source = source.replace(before, "if s[1] not in (2, 11):")
before = "table.append(row); self.symbols.append(row)"
assert source.count(before) == 1
source = source.replace(before, "table.append(row)\n                if s[1] == 2:\n                    self.symbols.append(row)")
exec(compile(source, str(basis), 'exec'))
