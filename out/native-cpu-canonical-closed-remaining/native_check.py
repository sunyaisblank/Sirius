"""Read-only native translation and full-function unwind-table inventory."""
from collections import Counter
from pathlib import Path
import bisect
import re
import subprocess
import sys

R = Path.cwd()
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf


def check(path, symbol):
    elf = Elf(path)
    functions = [row for row in elf.symbols if row['section'] and row['info'] & 15 == 2 and row['size']]
    assert len(functions) == 1
    function = functions[0]
    assert function['name'] == symbol and function['info'] == 18 and function['other'] & 3 == 2
    start, end = function['value'], function['value'] + function['size']
    assembly = subprocess.check_output(['/usr/bin/objdump', '-drw', str(path)], text=True)
    instructions = [(int(row[0], 16) - start, row[1]) for row in re.findall(
        r'^\s*([0-9a-f]+):\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9.]*)', assembly, re.M)]
    assert instructions and not [op for offset, op in instructions if op.startswith(('vfm', 'vfnm'))]
    relocs = Counter(re.findall(r'R_X86_64_PLT32\s+(sirius_local_value_\w+)-0x4', assembly))
    assert relocs
    frames = subprocess.check_output(['/usr/bin/readelf', '-W', '--debug-dump=frames-interp', str(path)], text=True)
    fdes = list(re.finditer(r'^([0-9a-f]+) .*\bFDE .*\bpc=([0-9a-f]+)\.\.([0-9a-f]+)\s*$', frames, re.M))
    matches = [i for i, row in enumerate(fdes) if (int(row[2], 16), int(row[3], 16)) == (start, end)]
    assert len(matches) == 1, 'One full-range FDE must cover the sole native entry.'
    i = matches[0]
    block = frames[fdes[i].end():fdes[i + 1].start() if i + 1 < len(fdes) else len(frames)]
    columns, rows = None, []
    for line in block.splitlines():
        if line.strip().startswith('LOC '):
            assert columns is None
            columns = line.split()
        elif columns and re.match(r'^[0-9a-f]+\s', line):
            fields = line.split()
            assert len(fields) == len(columns)
            row = dict(zip(columns, fields))
            offset = int(row.pop('LOC'), 16) - start
            assert 0 <= offset <= function['size']
            rows.append(dict(offset=offset, **row))
        elif rows and not line.strip():
            break
    assert rows and rows[0]['offset'] == 0
    assert all(a['offset'] < b['offset'] for a, b in zip(rows, rows[1:]))
    offsets = [row['offset'] for row in rows]
    def state(offset):
        i = bisect.bisect_right(offsets, offset) - 1
        assert i >= 0
        return rows[i]
    assert state(0)['CFA'] == 'rsp+8' and state(0)['ra'] == 'c-8'
    calls = [offset for offset, op in instructions if op in ('call', 'callq')]
    returns = [offset for offset, op in instructions if op in ('ret', 'retq')]
    assert calls and returns
    for offset in calls:
        row = state(offset)
        assert re.fullmatch(r'(?:rsp|rbp)\+[0-9]+', row['CFA']) and row['ra'] == 'c-8'
        if row['CFA'] == 'rbp+16':
            assert row['rbp'] == 'c-16'
    for offset in returns:
        assert state(offset)['CFA'] == 'rsp+8' and state(offset)['ra'] == 'c-8'
    return dict(instructions=len(instructions), emitted_fused_instructions=[], core_relocations=dict(relocs),
        sole_hidden_C_ABI_entry=symbol,
        unwind_table=dict(full_function_FDE=True, function_bytes=function['size'],
            relative_range=[0, function['size']], relative_CFI_rows=rows,
            covered_call_sites=len(calls), checked_return_offsets=returns,
            scope='Full FDE range, initial/return CFA and return-PC rule, and decoded CFA at every native call are checked. Saved-register rows are retained for final linked-code/prologue joins; actual same-thread job/parent/return unwind remains a required genuine C8 observation. This is not a runtime unwind or numerical verdict.'))
