"""Read-only exact code/relocation equality and reusable full-entry CFI check."""
from pathlib import Path
import bisect
import hashlib
import re
import subprocess
import sys

R = Path.cwd()
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf, pin

NAME = 'sirius_closed_dp_packet'
GP = dict(rbx='c-56', rbp='c-16', r12='c-48', r13='c-40', r14='c-32', r15='c-24')


def named(elf, symbol=NAME):
    rows = [s for s in elf.symbols if s['name'] == symbol and s['section'] and (s['info'] & 15) == 2 and s['size']]
    assert len(rows) == 1
    return rows[0]


def verify_unwind(path, symbol=NAME):
    """Function-relative CFI rows usable on the object and final linked DSO."""
    path = Path(path)
    elf = Elf(path)
    function = named(elf, symbol)
    start, end = function['value'], function['value'] + function['size']
    code = elf.function_bytes(function)
    assert code[0x1b:0x33] == bytes.fromhex('554889e54157415641554154534883e4c04881ec00c20300')
    frames = subprocess.check_output(['/usr/bin/readelf', '-W', '--debug-dump=frames-interp', str(path)], text=True)
    fdes = list(re.finditer(r'^([0-9a-f]+) .*\bFDE .*\bpc=([0-9a-f]+)\.\.([0-9a-f]+)\s*$', frames, re.M))
    matches = [i for i, match in enumerate(fdes) if (int(match[2], 16), int(match[3], 16)) == (start, end)]
    assert len(matches) == 1, 'One exact full-function FDE is required, not section presence.'
    index = matches[0]
    block = frames[fdes[index].end():fdes[index + 1].start() if index + 1 < len(fdes) else len(frames)]
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
    assert state(0x1b)['CFA'] == 'rsp+8'
    assert state(0x1c)['CFA'] == 'rsp+16' and state(0x1c)['rbp'] == 'c-16'
    assert state(0x1f)['CFA'] == 'rbp+16'
    assembly = subprocess.check_output(['/usr/bin/objdump', '-d', '--no-show-raw-insn',
        '--start-address=' + str(start), '--stop-address=' + str(end), str(path)], text=True)
    instructions = []
    for line in assembly.splitlines():
        match = re.match(r'^\s*([0-9a-f]+):\s+([a-z][a-z0-9.]*)\b(?:\s+(.*))?$', line)
        if match:
            instructions.append((int(match[1], 16) - start, match[2], match[3] or ''))
    calls = [offset for offset, op, _args in instructions if op in ('call', 'callq')]
    returns = [offset for offset, op, _args in instructions if op in ('ret', 'retq')]
    assert calls and 0x17bb in calls and returns
    for offset in calls:
        rules = state(offset)
        assert rules['CFA'] == 'rbp+16' and rules['ra'] == 'c-8'
        assert all(rules.get(reg) == rule for reg, rule in GP.items()), (offset, rules)
    for offset in returns:
        assert state(offset)['CFA'] == 'rsp+8' and state(offset)['ra'] == 'c-8', (offset, state(offset))
    assert all(row['ra'] == 'c-8' and row['CFA'] in ('rsp+8', 'rsp+16', 'rbp+16') for row in rows)
    return dict(full_function_FDE=True, symbol=symbol, function_bytes=function['size'],
        relative_range=[0, function['size']], relative_CFI_rows=rows,
        initial_CFA='rsp+8', body_CFA='rbp+16', complete_GP_save_rules=True,
        saved_GP_rules=GP, saved_return_PC='CFA-8', shrinkwrap_and_epilogues_checked=True,
        covered_call_sites=len(calls), checked_return_offsets=returns,
        same_thread_live_job_and_return_observation_still_required=True,
        scope='Decoded exact full-entry FDE/CFA/GP-save rules at all native calls and all returns, including pre-prologue zero-mask and restored epilogues. No runtime unwind/job/numerical success claim.')


def verify(original, candidate):
    original, candidate = Path(original), Path(candidate)
    a, b = Elf(original), Elf(candidate)
    assert a.kind == b.kind == 1 and '.eh_frame' not in a.names and '.debug_frame' not in a.names
    fa, fb = named(a), named(b)
    assert (fa['info'], fa['other'], fa['size'], fa['value']) == (fb['info'], fb['other'], fb['size'], fb['value'])
    assert fa['info'] == 18 and fa['other'] & 3 == 2
    for elf in (a, b):
        assert len([s for s in elf.symbols if s['section'] and (s['info'] & 15) == 2 and s['size']]) == 1
    assert a.function_bytes(fa) == b.function_bytes(fb)
    def relocs(elf):
        return sorted((elf.section_name(section), offset, kind, addend,
            target['name'], target['info'], target['other'], elf.section_name(target['section']), target['value'], target['size'])
            for section, offset, kind, addend, target in elf.relocations
            if elf.section_name(section) != '.eh_frame')
    assert relocs(a) == relocs(b)
    allocated = {}
    for i, section in enumerate(a.sections):
        if not section[2] & 2:
            continue
        j, = [j for j, name in enumerate(b.names) if name == a.names[i]]
        assert a.bytes(i) == b.bytes(j)
        allocated[a.names[i]] = dict(bytes=len(a.bytes(i)), sha256=hashlib.sha256(a.bytes(i)).hexdigest())
    assert {b.names[i] for i, s in enumerate(b.sections) if s[2] & 2} == set(allocated) | {'.eh_frame'}
    unwind = verify_unwind(candidate)
    return dict(passed=True, original=dict(path=str(original.relative_to(R)), **pin(original)),
        output=dict(path=str(candidate.relative_to(R)), **pin(candidate)),
        full_function_bytes_identical=True, normalized_code_relocations_identical=True,
        normalized_relocation_count=len(relocs(a)), all_original_allocated_sections_identical=allocated,
        sole_hidden_six_argument_C_ABI_preserved=True, unwind_table=unwind,
        numerical_qualification_transferred=False,
        scope='Only new unwind metadata is permitted. Function bytes, names/binding/visibility/size and every original non-unwind normalized relocation/allocated payload match exactly. Final DSO FDE must be joined separately before fresh genuine-job C8.')
