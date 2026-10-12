"""Gate exact artifact correspondence and join saved PCs to complete machine blocks."""
from pathlib import Path
from bisect import bisect_right
from collections import Counter, defaultdict
import hashlib
import json
import re
import resource
import subprocess
import sys
from elf_layout import elf

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
folder = root / 'bin/linux-gcc/portable-execution-attribution'
protected = root / 'attestations/software-vulkan/b5d2abf/jit-stage-linked-byte-binding/capture/jit-10.elf'
profile = root / 'attestations/software-vulkan/9624562/jit-dispatch-profile/corrected-pid-namespace'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
reference = elf(protected)
reports = []
for path in (folder / 'endpoint-machine-blocks.o', profile / 'jit-10.elf'):
    current = elf(path, '.Lsirius_machine_')
    assert current[2] == reference[2], (path, 'all function layouts')
    assert current[3] == reference[3], (path, 'all symbolic relocations')
    assert set(current[1]) == set(reference[1]), (path, 'all allocated sections')
    for name, before in reference[1].items():
        after = current[1][name]
        assert (before['flags'], before['align'], before['bytes'], before['data']) == (
            after['flags'], after['align'], after['bytes'], after['data']), (path, name)
    reports.append({'path': str(path.relative_to(root)), 'sha256': sha(path),
                    'complete_function_count': len(current[2]),
                    'complete_allocated_sections': len(current[1]),
                    'exact_symbolic_relocations': sum(map(len, current[3].values()))})

assembly = json.loads((work / 'machine-block-assembly.json').read_text())
assert sha(folder / 'endpoint.s') == assembly['assembly_sha256']
assert sha(folder / 'endpoint-machine-blocks.s') == assembly['annotated_assembly_sha256']
labels = elf(folder / 'endpoint-machine-blocks.o', '.Lsirius_machine_')[4]
blocks = assembly['blocks']
assert set(labels) == {b[key] for b in blocks for key in ('start_symbol', 'end_symbol')}
for block in blocks:
    section, start = labels[block['start_symbol']]
    end_section, end = labels[block['end_symbol']]
    function_section, function_start, size = reference[2][block['function']]
    assert section == end_section == function_section == '.ltext'
    assert function_start <= start <= end <= function_start + size
    block.update(section_start=start, section_end=end, function_start=start-function_start,
                 function_end=end-function_start)

disassembly = folder / 'endpoint-machine-disassembly.txt'
assert not disassembly.exists()
with disassembly.open('xb') as output:
    subprocess.run(['/usr/bin/objdump', '-d', '-w', str(folder / 'endpoint.o')],
                   stdout=output, check=True)
instructions = []
code = reference[1]['.ltext']['data']
next_address = 0
for line in disassembly.read_text().splitlines():
    match = re.match(r'^\s*([0-9a-f]+):\s*\t([0-9a-f ]+)\t(.+)$', line)
    if not match:
        continue
    address = int(match.group(1), 16)
    data = bytes.fromhex(match.group(2))
    operation = match.group(3).strip()
    assert address == next_address and data and data == code[address:address+len(data)]
    assert '(bad)' not in operation
    instructions.append({'section_offset': address, 'bytes': len(data), 'assembly': operation})
    next_address += len(data)
assert next_address == len(code), (next_address, len(code))
boundaries = {i['section_offset'] for i in instructions} | {len(code)}
for block in blocks:
    assert block['section_start'] in boundaries and block['section_end'] in boundaries, block
ordered = sorted((b for b in blocks if b['section_start'] < b['section_end']),
                 key=lambda b: b['section_start'])
for left, right in zip(ordered, ordered[1:]):
    assert left['section_end'] <= right['section_start']
starts = [b['section_start'] for b in ordered]
instruction_starts = [i['section_offset'] for i in instructions]
historical = elf(profile / 'jit-10.elf')
raw = json.loads((profile / 'jit-attribution.json').read_text())
object_record = next(o for o in raw['objects'] if o['file'] == 'jit-10.elf')
assert sha(profile / 'jit-10.elf') == object_record['sha256']
joined = []
weights = defaultdict(Counter)
dispositions = Counter()
for row in raw['mapped_instruction_samples']:
    if row['file'] != 'jit-10.elf':
        continue
    section, function_start, size = reference[2][row['function']]
    assert row['section_offset'] == function_start + row['function_offset']
    assert 0 <= row['function_offset'] < size
    assert row['ip'] == historical[1][section]['address'] + row['section_offset']
    index = bisect_right(starts, row['section_offset']) - 1
    block = ordered[index] if index >= 0 else None
    if block is not None and not (block['section_start'] <= row['section_offset'] < block['section_end']):
        block = None
    instruction = instructions[bisect_right(instruction_starts, row['section_offset']) - 1]
    assert instruction['section_offset'] <= row['section_offset'] < instruction['section_offset'] + instruction['bytes']
    status = 'machine_block' if block else 'alignment_or_interfunction_padding'
    dispositions[status] += row['count']
    aligned = instruction['section_offset'] == row['section_offset']
    dispositions['instruction_start' if aligned else 'instruction_interior'] += row['count']
    entry = {**row, 'disposition': status, 'instruction_start': aligned,
             'containing_instruction': instruction}
    if block:
        assert block['function'] == row['function']
        entry.update(machine_block=block['machine_block'], header_line=block['header_line'],
                     ir_comment=block['ir_comment'], block_start=block['section_start'],
                     block_end=block['section_end'])
        weights[(block['function'], block['machine_block'])][row['phase']] += row['count']
    joined.append(entry)
assert len(joined) == 383 and sum(row['count'] for row in joined) == 417
by_key = {(b['function'], b['machine_block']): b for b in blocks}
ranked = [{**by_key[key], 'weighted_samples': sum(counts.values()),
           'weighted_by_phase': dict(sorted(counts.items()))}
          for key, counts in sorted(weights.items(), key=lambda item: (-sum(item[1].values()), item[0]))]
scope = ('Exact historical Endpoint PC-offset to emitted machine-block correspondence only. '
         'All block headers and padding are accounted for. IR comments may reflect splitting, '
         'merging or motion; sampled PCs can have skid. No exclusive operation latency, '
         'current live-PC, unique producer-option, performance gain or release qualification claim.')
result = {'passed': True, 'scope': scope, 'reference_sha256': sha(protected),
          'artifact_gates': reports, 'assembly_sha256': assembly['assembly_sha256'],
          'objdump_sha256': sha(Path('/usr/bin/objdump')),
          'disassembly_sha256': sha(disassembly), 'whole_ltext_decoded_bytes': len(code),
          'whole_ltext_instruction_count': len(instructions), 'complete_emitted_blocks': len(blocks),
          'retained_headers': sum(b['header_kind']=='retained_label' for b in blocks),
          'fallthrough_headers': sum(b['header_kind']=='fallthrough_comment' for b in blocks),
          'empty_blocks': sum(b['section_start']==b['section_end'] for b in blocks),
          'historical_sample_rows': len(joined), 'historical_weighted_samples': 417,
          'sample_dispositions': dict(dispositions), 'sampled_machine_blocks': len(ranked),
          'top_blocks': ranked[:20]}
for name, value in [('machine-block-offsets.json', blocks), ('historical-machine-block-samples.json', joined),
                    ('historical-machine-block-ranking.json', ranked), ('machine-block-verification.json', result)]:
    path = work / name
    assert not path.exists()
    path.write_text(json.dumps(value, indent=2) + '\n')
print(json.dumps({k:v for k,v in result.items() if k != 'top_blocks'}, indent=2))
