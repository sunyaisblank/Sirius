"""Join pinned private inputs, unwind metadata and callback PCs to the final ELF."""

from collections import Counter
from pathlib import Path
import re
import struct
import subprocess
import sys

from prepare_provider import AX, R
from resume_raw import pin, read_json, write_json
from run_native_context import linked_symbols
from run_value_abi import expect_pin

sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf, associate


STAGES = {'transport': 'transport', 'endpoint': 'endpoint', 'dopri_phase': 'dp'}


def named(elf, symbol):
    row, = (row for row in elf.symbols if row['name'] == symbol and row['section']
            and row['size'] and row['info'] & 15 == 2)
    return row


def disassemble(elf, name):
    symbol = named(elf, name)
    lines = subprocess.check_output(['/usr/bin/objdump', '-d', '-Mintel', f'--disassemble={name}',
                                     str(elf.path)], text=True).splitlines()
    return symbol['value'], [(int(match[1], 16), match[2].strip()) for line in lines
                            if (match := re.match(r'^\s*([\da-f]+):\s+(?:[\da-f]{2}\s+)+\s*(\S.*)', line))]


def call_return(instructions, target):
    return_pc, = (after for (_, op), (after, _) in zip(instructions, instructions[1:])
                  if re.fullmatch(r'call\s+[\da-f]+ <' + re.escape(target) + '>', op))
    return return_pc


def fde_rows(table, address, size):
    matches = list(re.finditer(r'^[0-9a-f]+ .*\bFDE .*\bpc=([0-9a-f]+)\.\.([0-9a-f]+)\s*$', table, re.M))
    index, = (index for index, match in enumerate(matches)
              if (int(match[1], 16), int(match[2], 16)) == (address, address + size))
    next_start = matches[index + 1].start() if index + 1 < len(matches) else len(table)
    columns, rows = None, []
    for line in table[matches[index].end():next_start].splitlines():
        if line.strip().startswith('LOC '):
            columns = line.split()
        elif columns and re.match(r'^[0-9a-f]+\s', line):
            row = dict(zip(columns, line.split(), strict=True))
            rows.append(dict(offset=int(row.pop('LOC'), 16) - address, **row))
        elif rows and not line.strip():
            break
    return rows


def main():
    scratch, provider = (Path(arg).resolve() for arg in sys.argv[1:])
    component = read_json(scratch / 'provider-component.json')
    native = read_json(scratch / 'native-replay.json')['results']
    output = R / component['output']['path']
    expect_pin(output, component['output'])
    expect_pin(R / native['link']['output']['path'], native['link']['output'])
    local_core_names = set(linked_symbols(R / native['link']['output']['path'])['defined_symbols'])
    local_core_names = {name for name in local_core_names if name.startswith('sirius_local_value_')}

    # Only declared static relocation operands may differ from the pinned objects.
    association = associate(provider, [*read_json(AX)['native_objects'], component['output']])
    elf = Elf(provider)
    elf.path = provider
    wide = Elf(output)
    relocations = Counter()
    wide_prefixes = ('sirius_wide_strict_', 'sirius_wide_value_', 'sirius_local_value_')
    for function in association['objects'][-1]['functions']:
        old = named(wide, function['symbol'])
        new = named(elf, function['symbol'])
        code = elf.function_bytes(new)
        if new['info'] >> 4 != 0:
            raise RuntimeError(f'Retained native symbol is exported: {function["symbol"]}')
        for section, location, kind, addend, target in wide.relocations:
            if section != old['section'] or not old['value'] <= location < old['value'] + old['size']:
                continue
            if not target['name'].startswith(wide_prefixes):
                continue
            offset = location - old['value']
            if kind != 4 or addend != -4 or code[offset - 1] != 0xE8:
                raise RuntimeError('Expected a direct native call relocation')
            destination = new['value'] + offset + 4 + struct.unpack_from('<i', code, offset)[0]
            if destination != named(elf, target['name'])['value']:
                raise RuntimeError(f'Native call does not reach {target["name"]}')
            relocations[target['name']] += 1
    assert sum(relocations.values()) == 8247 and len(set(relocations) & local_core_names) == 9

    table = subprocess.check_output(['/usr/bin/readelf', '-W', '--debug-dump=frames-interp', str(provider)], text=True)
    unwind = {}
    for stage, slug in STAGES.items():
        symbol = named(elf, f'sirius_wide_strict_{slug}_packet')
        rows = fde_rows(table, symbol['value'], symbol['size'])
        original = native['emit-' + stage]['native']['unwind_table']
        assert rows == original['relative_CFI_rows']
        unwind[stage] = dict(address=symbol['value'], bytes=symbol['size'], cfi_rows=len(rows),
                             relative_cfi_matches_object=True)

    callback_name, mode_name = 'cs_exec_fn', 'lp_native_cpu_packet_mode'
    _, callback = disassemble(elf, callback_name)
    _, mode = disassemble(elf, mode_name)
    mode_return = call_return(callback, mode_name)
    before = [op for address, op in callback if address < mode_return]
    assert [op for op in before if re.match(r'\w+\s+(?:rbx|ebx|bx|bl),', op)] == ['mov    rbx,rdi']
    assert before[-5:-1] == ['mov    ecx,r14d', 'mov    edx,r13d', 'mov    esi,ebp', 'mov    rdi,rbx']
    assert before[-6] == 'xor    r8d,r8d'
    assert 'mov    DWORD PTR [rsp+0x4],eax' in before
    assert 'and    edi,0xffff003f' in before and 'or     edi,0x1f80' in before
    assert [op for address, op in callback if mode_return <= address < mode_return + 24][-2:] == [
        'mov    edi,DWORD PTR [rsp+0x4]', next(op for address, op in callback
        if address > mode_return and re.search(r'call\s+[\da-f]+ <util_fpstate_set>', op))]
    assert any('stmxcsr' in op for _, op in disassemble(elf, 'util_fpstate_get')[1])
    assert any('ldmxcsr' in op for _, op in disassemble(elf, 'util_fpstate_set')[1])
    packets = {}
    for stage in STAGES:
        for prefix in ('', 'wide_'):
            name = f'sirius_native_{stage}_rows_{prefix}admitted_packet'
            packets[name] = {'entry': named(elf, name)['value'], 'return': call_return(mode, name)}

    write_json(scratch / 'linked-provider.json', {
        'provider': {'path': str(provider.relative_to(R)), **pin(provider)},
        'native_opcode_associations': len(association['objects']),
        'direct_wide_join_relocations': dict(sorted(relocations.items())), 'unwind': unwind,
        'callback': {'mode_return': mode_return, 'job_register': 'rbx', 'saved_mxcsr_offset': 4},
        'packet_entries': packets, 'provider_component': component['output'],
    })


if __name__ == '__main__':
    main()
