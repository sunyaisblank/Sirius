"""Lower every admitted canonical layer using unchanged generated arithmetic.

Default is read-only. --write creates missing artifacts, never replaces them.
No shader generation, compilation, provider mutation, or execution occurs here.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
basis = root / 'bin/linux-gcc/native-cpu-owner-context/typed-batch'
tables = root / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch'
destination = root / 'bin/linux-gcc/native-cpu-static-program/typed-batch'
assert sys.argv[1:] in ([], ['--write'])
write = sys.argv[1:] == ['--write']


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert write, 'Missing derived input: ' + str(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


def function(text, name):
    start = text.index('static bool ' + name + '(')
    opening = text.index('{', start)
    depth = 1
    end = opening + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return start, opening, end


specs = {
    'transport': {'layers': [117, 42],
                  'args': 'programBase_1, count_0, registers_1, layer_0, start_0, end_0, row_3, rowBase_7, kernelContext_12'},
    'endpoint': {'layers': [111, 84],
                 'args': 'programBase_2, count_1, registers_2, layer_0, start_0, end_0, row_1, rowBase_10, kernelContext_12'},
    'dopri_phase': {'layers': [45],
                    'args': 'programBase_1, count_0, registers_1, layer_0, start_0, end_0, row_2, rowBase_6, kernelContext_8'},
}
records = {}
total_instructions = total_layers = 0
for slug, spec in specs.items():
    source = basis / (slug + '.native.hpp')
    original = source.read_text()
    table = tables / (slug + '.program.hpp')
    table_text = table.read_text()
    adapter = tables / (slug + '.cpp')
    stride_match = re.search(r'first_x, y, z, grid_x, grid_y, grid_z, groups, (\d+), (\d+),', adapter.read_text())
    assert stride_match
    stride = int(stride_match[1])
    values = [int(v) for v in re.findall(r'(\d+)u', table_text)]
    start, opening, end = function(original, 'EvaluateOne_0')
    signature = original[start:opening].strip()
    params = signature[signature.index('(') + 1:-1].split(', ')
    assert len(params) == 6
    index, row, rowbase, programbase, registers, context = [p.split()[-1] for p in params]
    body = original[opening:end]
    bindings = re.findall(r'uint32_t (\w+) = [^;]*inputs_0\.Load\(address_0(?: \+ [1-4]U)?\);', body)
    assert bindings == ['op_0', 'destination_0', 'a_31', 'b_22' if slug != 'endpoint' else 'b_23', 'c_4' if slug == 'endpoint' else 'c_2'], bindings
    last = re.search(r'    uint32_t ' + bindings[-1] + r' = [^;]+;', body)
    tail = body[last.end():]
    assert not any(re.search(r'\b' + re.escape(p) + r'\b', tail) for p in (index, programbase, 'address_0'))
    static = ('template<uint32_t StaticOp>\nstatic bool EvaluateOneStatic_0(' +
              ', '.join('uint32_t ' + p for p in bindings[1:]) + ', ' +
              ', '.join(params[i] for i in (1, 2, 4, 5)) + ')\n{\n' +
              '    static_assert(StaticOp <= 11U);\n    constexpr uint32_t op_0 = StaticOp;' + tail + '\n')
    assert static[static.index(tail):static.index(tail) + len(tail)] == tail
    helper = '\n// Canonical instruction metadata is admitted before any native writes.\n' + static
    dispatches = []
    program_records = []
    cursor = 0
    for family, layer_count in enumerate(spec['layers']):
        base = cursor
        count, reg_count, roots = values[base:base + 3]
        instruction_base = base + 3 + roots
        offsets_base = instruction_base + 5 * count
        offsets = values[offsets_base:offsets_base + layer_count + 1]
        assert len(offsets) == layer_count + 1 and offsets[0] == 0 and offsets[-1] == count
        assert all(0 < b - a <= 64 for a, b in zip(offsets, offsets[1:]))
        tuples = [values[instruction_base + 5*i:instruction_base + 5*i + 5] for i in range(count)]
        assert all(len(t) == 5 and t[0] <= 11 and t[1] < reg_count for t in tuples)
        helper += f'static_assert(sirius_cpu_bridge::{slug}_program[{base}] == {count}U);\n'
        cases = []
        emitted = []
        for layer, (first, limit) in enumerate(zip(offsets, offsets[1:])):
            name = f'EvaluateStaticLayer_{family}_{layer}_0'
            helper += (f'static void {name}(uint32_t {row}, uint32_t {rowbase}, uint32_t {registers}, KernelContext_0 * {context})\n{{\n')
            for instruction in range(first, limit):
                op, dst, a, b, c = tuples[instruction]
                # All dependencies execute in original index order, even after a fault.
                helper += (f'    if(!EvaluateOneStatic_0<{op}U>({dst}U, {a}U, {b}U, {c}U, {row}, {rowbase}, {registers}, {context}))\n'
                           f'        RejectRetainedProgram_0({context});\n')
                emitted.append(tuples[instruction])
            helper += '}\n'
            cases.append(f'        case {layer}U:\n            if(start != {first}U || end != {limit}U) return false;\n            {name}(row, rowBase, registers, context);\n            return true;\n')
        assert emitted == tuples
        dispatches.append(f'    if(programBase == tableBase + {base}U && count == {count}U && registers == {reg_count}U)\n    {{\n        switch(layer)\n        {{\n' + ''.join(cases) + '        default: return false;\n        }\n    }\n')
        cursor = offsets_base + layer_count + 1
        program_records.append({'base': base, 'instructions': count, 'registers': reg_count,
                                'roots': roots, 'layers': layer_count,
                                'canonical_tuples': pin(json.dumps(tuples, separators=(',', ':')).encode()),
                                'opcodes': dict(sorted(Counter(t[0] for t in tuples).items())),
                                'all_original_tuples_in_original_layer_order': True})
        total_instructions += count
        total_layers += layer_count
    assert cursor == len(values), (slug, cursor, len(values))
    helper += ('static bool EvaluateKnownLayer_0(uint32_t programBase, uint32_t count, uint32_t registers,\n'
               '                                 uint32_t layer, uint32_t start, uint32_t end,\n'
               '                                 uint32_t row, uint32_t rowBase, KernelContext_0 * context)\n{\n'
               f'    const uint32_t tableBase = 1U + context->globalParams_0->inputs_0.Load(0U) * {stride}U;\n' +
               ''.join(dispatches) + '    return false;\n}\n\n')
    # The existing i=start assignment, outer layer checks and publication stay intact.
    program_start, _, program_end = function(original, 'EvaluateProgram_0')
    program_body = original[program_start:program_end]
    call = re.search(r'^ +bool \w+ = EvaluateOne_0\([^;]+;', program_body, re.M)
    assert call
    loop_start = program_body.rfind('            for(;;)', 0, call.start())
    # Endpoint nests the same loop one indentation level deeper.
    loop_start = program_body.rfind('\n', 0, loop_start) + 1
    brace = program_body.index('{', loop_start)
    depth, loop_end = 1, brace + 1
    while depth:
        depth += (program_body[loop_end] == '{') - (program_body[loop_end] == '}')
        loop_end += 1
    old_loop = program_body[loop_start:loop_end]
    assert old_loop.count('EvaluateOne_0(') == 1 and call.group().strip() in old_loop
    indent = re.match(r' *', old_loop).group()
    wrapped = (indent + 'if(!EvaluateKnownLayer_0(' + spec['args'] + '))\n' + indent + '{\n' +
               ''.join('    ' + line if line.strip() else line for line in old_loop.splitlines(True)) +
               '\n' + indent + '}')
    changed_program = program_body[:loop_start] + wrapped + program_body[loop_end:]
    assert changed_program.replace(wrapped, old_loop, 1) == program_body
    derived = (f'#include <{slug}.program.hpp>\n' + original[:program_start] + helper +
               changed_program + original[program_end:])
    assert original[start:end] in derived
    path = destination / (slug + '.native.hpp')
    publish(path, derived.encode())
    records[slug] = {'basis_path': str(source.relative_to(root)), 'basis': pin(source.read_bytes()),
                     'table_path': str(table.relative_to(root)), 'table': pin(table.read_bytes()),
                     'adapter_path': str(adapter.relative_to(root)), 'adapter': pin(adapter.read_bytes()),
                     'input_row_words_from_original_packet_guard': stride,
                     'derived_path': str(path.relative_to(root)), 'derived': pin(derived.encode()),
                     'programs': program_records, 'original_arithmetic_tail': pin(tail.encode()),
                     'arithmetic_tail_byte_identical': True, 'generic_fallback_byte_identical': True,
                     'outer_layer_validation_and_publication_byte_identical': True,
                     'other_generated_header_bytes_changed': False}
assert total_instructions == 9380 and total_layers == 399
record = {'scope': 'Complete canonical native lowering only. At most12 opcode template bodies per stage; operands are literal uint32 arguments, no instruction metadata reads in static helpers. Original generated arithmetic tail, generic fallback, outer layer checks, root publication, context lifetime and64logicalGroupcalls retained. All instructions execute in original layer/index order, including after structural false; invalid numerical results still write all5words. No compile/runtime or benefit claim.',
          'instructions': total_instructions, 'layers': total_layers, 'headers': records}
publish(work / 'header-preparation.json', (json.dumps(record, indent=2) + '\n').encode())
print(json.dumps({'derived_headers_checked': 3, 'instructions': total_instructions,
                  'layers': total_layers, 'arithmetic_changed': False, 'compiled_or_executed': False}))
