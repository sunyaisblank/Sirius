"""Pack independent canonical instructions; preserve every original layer and value flow.

This emits data for a future native field-SIMD executor. It does not compile or
execute math, change the canonical wire program, or qualify any runtime.
"""
from collections import Counter, defaultdict
from pathlib import Path
import hashlib
import json
import re

R = Path.cwd()
W = Path(__file__).resolve().parent
DEST = R / 'bin/linux-gcc/native-cpu-field-simd/generated'
ARITY = {0: 0, 1: 0, 2: 2, 3: 2, 4: 2, 5: 2, 6: 1, 7: 1, 8: 1, 9: 2, 10: 3, 11: 3}


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write(data)


def helper(text):
    found = re.search(r'static bool EvaluateOneStatic_0\(uniform uint32_t (\w+), uniform uint32_t (\w+), uniform uint32_t (\w+), uniform uint32_t (\w+),', text)
    assert found
    opening = text.index('{', found.end()); end = opening + 1; depth = 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}'); end += 1
    body = text[opening:end]
    dst, a, b, c = found.groups()
    op, = re.findall(r'uniform const uint32_t (\w+) = StaticOp;', body)
    reads = re.findall(r'ReadValueUniform_0\((\w+),', body)
    assert reads == [a, b, c]
    positions = {n: body.index(f'{op} == {n}U') for n in range(12)}
    ra, rb, rc = [body.index('ReadValueUniform_0(' + name + ',') for name in (a, b, c)]
    assert positions[0] < positions[1] < ra < positions[6] < positions[7] < positions[8] < rb
    assert rb < positions[2] < positions[3] < positions[4] < positions[5] < positions[9] < positions[10] < positions[11] < rc
    assert body.count('WriteValueUniform_0(') == 1 and 'WriteValueUniform_0(' + dst + ',' in body
    segment = body[positions[1]:ra]
    limit, = re.findall(r'if\(' + re.escape(a) + r' >= (\d+)U\)', segment)
    return {'source': pin(body.encode()), 'input_limit': int(limit),
            'scratch_read_roles_in_source_order': reads, 'destination_role': dst,
            'source_branch_shape_checked_for_arity': True}


def reads(instruction):
    op, _dst, a, b, c = instruction
    return (a, b, c)[:ARITY[op]]


def dependencies(instructions, first, end):
    # Include RAW, WAR and WAW, since canonical physical registers are reused.
    writers, readers, result = {}, defaultdict(set), {}
    for index in range(first, end):
        instruction = instructions[index]; dst = instruction[1]; prior = set()
        for register in reads(instruction):
            if register in writers:
                prior.add(writers[register])
            readers[register].add(index)
        if dst in writers:
            prior.add(writers[dst])
        prior.update(readers[dst] - {index})
        result[index] = prior
        writers[dst] = index; readers[dst].clear()
    return result


def schedule(instructions, first, end, width):
    if width == 1:
        return [[index] for index in range(first, end)]
    pending = dependencies(instructions, first, end)
    done, groups = set(), []
    while pending:
        available = defaultdict(list)
        for index, prior in pending.items():
            if prior <= done:
                available[instructions[index][0]].append(index)
        assert available, 'cycle or missing dependency'
        # Deterministic opcode-coherent packs; no instruction or row is dropped.
        op = min(available, key=lambda key: (-len(available[key]), min(available[key]), key))
        group = sorted(available[op])[:width]
        assert all(instructions[index][0] == op for index in group)
        groups.append(group)
        for index in group:
            del pending[index]
        done.update(group)
    return groups


def validate(instructions, first, end, groups, registers):
    flat = [index for group in groups for index in group]
    assert sorted(flat) == list(range(first, end)) and len(flat) == len(set(flat))
    # Independently replay symbolic register definitions in serial and packed
    # orders. Check each node's ordered operands and every final scratch cell.
    initial = [('entry', register) for register in range(registers)]
    serial, packed = list(initial), list(initial)
    expected = {}
    for index in range(first, end):
        expected[index] = tuple(serial[register] for register in reads(instructions[index]))
        serial[instructions[index][1]] = ('node', index)
    for group in groups:
        destinations = [instructions[index][1] for index in group]
        assert len(destinations) == len(set(destinations))
        assert all(not (set(reads(instructions[index])) & (set(destinations) - {instructions[index][1]})) for index in group)
        observed = {index: tuple(packed[register] for register in reads(instructions[index])) for index in group}
        assert all(observed[index] == expected[index] for index in group)
        for index in group:
            packed[instructions[index][1]] = ('node', index)
    assert packed == serial


def array(name, values):
    rows = [', '.join(str(value) + 'U' for value in values[start:start + 16]) for start in range(0, len(values), 16)]
    return 'static const uniform uint32_t ' + name + '[' + str(len(values)) + '] = {\n    ' + ',\n    '.join(rows) + '\n};\n'


source_manifest = R / 'out/native-cpu-static-program/header-preparation.json'
metadata = json.loads(source_manifest.read_text())
record = {'status': 'static schedule compiler preparation only; no native object or math workload',
          'source_inputs': {str(source_manifest.relative_to(R)): pin(source_manifest.read_bytes())},
          'families': {}, 'instructions': 0, 'original_layers': 0,
          'scope': 'Preserve all canonical tuples, original layer boundaries, ordered operand values and every final scratch cell. Packing changes only independent instruction order within a layer, with explicit RAW/WAR/WAW dependencies and an independent symbolic replay. Same arithmetic bodies, original roots/publication/refusal and runtime address/FP ownership still require implementation and full original acceptance; this is not runtime or scientific qualification.'}
for slug in ('transport', 'endpoint', 'dopri_phase'):
    table = R / metadata['headers'][slug]['table_path']
    source = R / ('bin/linux-gcc/native-cpu-zero-term/ispc-source/' + slug + '.ispc')
    scratch_words = int(re.search(r'retainedScratch_0\[(\d+)\]', source.read_text())[1])
    words = [int(value) for value in re.findall(r'(\d+)u', table.read_text())]
    shape = helper(source.read_text())
    assert pin(table.read_bytes()) == metadata['headers'][slug]['table']
    record['source_inputs'][str(table.relative_to(R))] = pin(table.read_bytes())
    record['source_inputs'][str(source.relative_to(R))] = pin(source.read_bytes())
    stage, emitted, cursor = [], '// Canonical native field schedules; arithmetic is unchanged.\n', 0
    for family, old in enumerate(metadata['headers'][slug]['programs']):
        assert old['base'] == cursor
        count, registers, roots = words[cursor:cursor + 3]
        assert (count, registers, roots) == (old['instructions'], old['registers'], old['roots'])
        assert registers * 5 == scratch_words
        begin = cursor + 3 + roots
        instructions = [tuple(words[begin + 5*i:begin + 5*i + 5]) for i in range(count)]
        offsets = words[begin + 5*count:begin + 5*count + old['layers'] + 1]
        assert offsets[0] == 0 and offsets[-1] == count and all(0 < b-a <= 64 for a,b in zip(offsets, offsets[1:]))
        assert all(op in ARITY and dst < registers and all(src < registers for src in reads(row)) and
                   (op != 1 or a < shape['input_limit']) for row in instructions for op,dst,a,_b,_c in [row])
        assert pin(json.dumps([list(row) for row in instructions], separators=(',', ':')).encode()) == old['canonical_tuples']
        family_record = {'base': cursor, 'instructions': count, 'registers': registers,
                         'roots': roots, 'layers': old['layers'], 'source_arity_check': shape, 'tiles': {}}
        for width in (1, 2, 4, 8):
            groups, layer_packs = [], [0]
            for first, end in zip(offsets, offsets[1:]):
                packed = schedule(instructions, first, end, width)
                validate(instructions, first, end, packed, registers)
                groups.extend(packed); layer_packs.append(len(groups))
            assert sum(map(len, groups)) == count
            prefix = f'SiriusField_{slug}_{family}_{width}'
            emitted += array(prefix + '_LayerPacks', layer_packs)
            emitted += array(prefix + '_Opcode', [instructions[group[0]][0] for group in groups])
            emitted += array(prefix + '_ActiveNodes', list(map(len, groups)))
            emitted += array(prefix + '_Instruction', [index for group in groups for index in group + [0xffffffff]*(width-len(group))])
            family_record['tiles'][str(width)] = {'packs': len(groups), 'active_nodes': count,
                'field_lane_utilization': count / (len(groups)*width),
                'pack_size_histogram': dict(sorted(Counter(map(len, groups)).items())),
                'all_original_instructions_once': True, 'ordered_value_flow_and_final_scratch_equal': True}
        stage.append(family_record); record['instructions'] += count; record['original_layers'] += old['layers']
        cursor = begin + 5*count + old['layers'] + 1
    assert cursor == len(words)
    path = DEST / (slug + '.ispc.inc')
    publish(path, emitted.encode())
    record['families'][slug] = {'programs': stage, 'generated_path': str(path.relative_to(R)), 'generated': pin(emitted.encode())}
assert record['instructions'] == 9380 and record['original_layers'] == 399
publish(W / 'schedule-preparation.json', (json.dumps(record, indent=2) + '\n').encode())
print(json.dumps({'generated_stages': 3, 'original_instructions': 9380, 'original_layers': 399,
                  'field_widths': [1, 2, 4, 8], 'independent_symbolic_replay': 'PASS', 'workloads_executed': 0}))
