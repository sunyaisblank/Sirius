"""Join ordered live pure expressions to original serialized program writes.

This reads the unchanged producer and complete generated tables. It captures
allocator locals at return without changing the producer or allocation plan.
It proves a source/table mapping, not numeric equality or performance.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib
import json
import re
import struct

ROOT = Path('/home/astra/.project/Sirius')
WORK = Path(__file__).resolve().parent
SOURCE = ROOT / 'src/sirius/kernels/retained_program.py'
HEADER = ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
SOURCE_SHA = '3754f0911b2b1a3d58276b23f4468d61a6092c65a1b522d7b557e373c6e96516'
HEADER_SHA = '3f7aca3372bb872934e978fb0270f2375ce8d2361555b35c8b891e42e07161cc'
HEAD = '9fbd0e8b0d549a13476d066f2877dcd2cf9a94d9'
PURE_INPUTS = frozenset([0, 1, 2, 3, 4, 5, 6, 7, 44])

def seal(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def dependencies(operation):
    op, a, b, c = operation
    assert 0 <= op <= 11
    return () if op in (0, 1) else ((a,) + (() if op in (6, 7, 8) else (b,))
                                   + ((c,) if op in (10, 11) else ()))

def closure(operations, roots):
    seen, todo = set(), list(roots)
    while todo:
        index = todo.pop()
        if index not in seen:
            seen.add(index)
            todo.extend(dependencies(operations[index]))
    return seen

source, header = SOURCE.read_bytes(), HEADER.read_bytes()
assert seal(source)['sha256'] == SOURCE_SHA
assert seal(header)['sha256'] == HEADER_SHA
namespace = {'__name__': 'readonly_original_program_producer'}
exec(compile(source, str(SOURCE), 'exec'), namespace)
assert sys.getprofile() is None
original_compile = namespace['compile_program']
allocator_code = namespace['compile_parallel_program'].__code__

def build(name):
    captured = {}
    def capture_roots(outputs, *args, **kwargs):
        captured['roots'] = tuple(outputs)
        return original_compile(outputs, *args, **kwargs)
    def observe(frame, event, argument):
        if frame.f_code is allocator_code and event == 'return':
            assert 'layers' not in captured
            captured['layers'] = tuple(tuple(layer) for layer in frame.f_locals['layers'])
            captured['slots'] = dict(frame.f_locals['slots'])
            captured['live'] = frozenset(frame.f_locals['live'])
    namespace['compile_program'] = capture_roots
    try:
        sys.setprofile(observe)
        plan = namespace[name](True)
    finally:
        sys.setprofile(None)
        namespace['compile_program'] = original_compile
    operations = tuple(namespace['ops'])
    assert captured['live'] == closure(operations, captured['roots'])
    ordered = [index for layer in captured['layers']
               for index in sorted(layer, key=lambda i: (operations[i][0], i))]
    assert len(ordered) == plan['instructions'] == len(set(ordered))
    assert set(ordered) == captured['live']
    assert [0] + [sum(len(layer) for layer in captured['layers'][:end])
                  for end in range(1, len(captured['layers']) + 1)] == plan['layer_offsets']
    assert [captured['slots'][i] for i in captured['roots']] == plan['outputs']
    writes = {}
    layer_of = {index: layer for layer, nodes in enumerate(captured['layers']) for index in nodes}
    for ordinal, index in enumerate(ordered):
        op, a, b, c = operations[index]
        slots = captured['slots']
        if op not in (0, 1):
            a, b, c = slots[a], (0 if op in (6, 7, 8) else slots[b]), (slots[c] if op in (10, 11) else 0)
        expected = (op, slots[index], a, b, c)
        actual = tuple(plan['operations'][5 * ordinal:5 * ordinal + 5])
        assert actual == expected, (name, index, ordinal, actual, expected)
        writes[index] = {'source_node': index, 'instruction_ordinal': ordinal,
                         'layer': layer_of[index], 'tuple': list(actual), 'destination': actual[1]}
    owners = {}
    for index in reversed(ordered):
        destination = writes[index]['destination']
        following = owners.get(destination)
        writes[index]['next_write_ordinal'] = None if following is None else writes[following]['instruction_ordinal']
        writes[index]['final_register_owner'] = following is None
        owners[destination] = index
    return {'operations': operations, 'roots': captured['roots'], 'plan': plan,
            'live': captured['live'], 'writes': writes}

graphs = {}
for family, transport, endpoint in [
    ('general', 'build_transport_program', 'build_endpoint_program'),
    ('Schwarzschild', 'build_schwarzschild_transport_program', 'build_schwarzschild_endpoint_program')]:
    graphs[family] = {'Transport': build(transport), 'Endpoint': build(endpoint)}

arrays = {}
pattern = r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};'
for count, name, body in re.findall(pattern, header.decode(), re.S):
    values = tuple(int(word[:-1], 0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u', body))
    assert len(values) == int(count) and name not in arrays
    arrays[name] = values
assert len(arrays) == 40

def encode(plan):
    return tuple([plan['instructions'], plan['registers'], len(plan['outputs'])]
                 + plan['outputs'] + plan['operations'] + plan['layer_offsets'])

table_records = {}
for stage, header_words, packet_stride in [('Transport', 43, 230), ('Endpoint', 103, 225)]:
    combined = ()
    for family in ['general', 'Schwarzschild']:
        graph = graphs[family][stage]
        encoded = encode(graph['plan'])
        graph['table_offset'] = len(combined)
        graph['header_words'] = header_words
        graph['packet_stride'] = packet_stride
        graph['encoded'] = encoded
        combined += encoded
    assert combined == arrays['k' + stage + 'Program'], stage
    data = struct.pack('<' + 'I' * len(combined), *combined)
    table_records[stage] = {'array': 'k' + stage + 'Program', **seal(data),
                            'complete_original_table_exact': True, 'packet_row_stride_words': packet_stride,
                            'shader_program_address': '1 + packet_row_stride_words * stored_rows + selected_family_offset',
                            'instruction_address_offset_words': header_words}

results = {}
for family, stages in graphs.items():
    intern, maps = {}, {}
    for stage, graph in stages.items():
        ids, input_sets = [], []
        for index, operation in enumerate(graph['operations']):
            op, a, b, c = operation
            deps = dependencies(operation)
            assert all(dependency < index for dependency in deps)
            key = (op, a) if op in (0, 1) else (op, tuple(ids[dependency] for dependency in deps))
            if key not in intern:
                intern[key] = len(intern)
            ids.append(intern[key])
            inputs = (frozenset([a]) if op == 1 else frozenset()) if op in (0, 1) else frozenset().union(*(input_sets[d] for d in deps))
            input_sets.append(inputs)
        eligible = sorted(index for index in graph['live'] if input_sets[index] <= PURE_INPUTS)
        matching = {ids[index]: index for index in eligible}
        assert len(matching) == len(eligible), ('duplicate structural node', family, stage)
        maps[stage] = matching
        graph['input_sets'] = input_sets
        graph['structural_ids'] = ids
    common = maps['Transport'].keys() & maps['Endpoint'].keys()
    pairs = sorted((maps['Transport'][key], maps['Endpoint'][key]) for key in common)
    records = []
    for transport, endpoint in pairs:
        record = {'opcode': stages['Transport']['operations'][transport][0],
                  'input_records': sorted(stages['Transport']['input_sets'][transport])}
        assert record['input_records'] == sorted(stages['Endpoint']['input_sets'][endpoint])
        for stage, index in [('Transport', transport), ('Endpoint', endpoint)]:
            graph = stages[stage]
            write = dict(graph['writes'][index])
            offset = graph['table_offset'] + graph['header_words'] + 5 * write['instruction_ordinal']
            assert list(arrays['k' + stage + 'Program'][offset:offset + 5]) == write['tuple']
            write['combined_table_tuple_offset_words'] = offset
            write['scratch_component_indices'] = [5 * write['destination'] + i for i in range(5)]
            if stage == 'Endpoint':
                write['in_original_metric_tangent_prefix'] = write['instruction_ordinal'] < graph['plan']['prefix_instructions']
            record[stage] = write
        records.append(record)
    arithmetic = [record for record in records if record['opcode'] not in (0, 1)]
    # This is the address/write mapping of identical source expressions; actual
    # supplied tables and all normalized five-word inputs remain runtime keys.
    results[family] = {
        'common_live_nodes': len(records), 'common_arithmetic_nodes': len(arithmetic),
        'opcode_counts': {str(op): sum(record['opcode'] == op for record in records) for op in range(12)},
        'arithmetic_transport_registers_recycled_later': sum(not record['Transport']['final_register_owner'] for record in arithmetic),
        'arithmetic_endpoint_prefix_nodes': sum(record['Endpoint']['in_original_metric_tangent_prefix'] for record in arithmetic),
        'programs': {stage: {'live_nodes': len(graph['live']), 'registers': graph['plan']['registers'],
                            'layers': len(graph['plan']['layer_offsets']) - 1,
                            'selected_family_table_offset_words': graph['table_offset'],
                            'original_tuple_mappings_verified': len(graph['writes']),
                            'prefix_instructions': graph['plan'].get('prefix_instructions')}
                     for stage, graph in stages.items()},
        'matches': records}

assert SOURCE.read_bytes() == source and HEADER.read_bytes() == header
destination = WORK / 'mapping-result.json'
assert not destination.exists()
destination.write_text(json.dumps({'accepted_source_mapping': True, 'source_revision': HEAD,
    'source': seal(source), 'header': seal(header), 'complete_tables': table_records, 'families': results,
    'normalized_input_key_records': sorted(PURE_INPUTS),
    'runtime_equality_proved': False, 'owner_interface_implemented': False, 'performance_trial_authorized': False,
    'scope': 'Exact ordered source-node to complete original table/write mapping only. No source algebra, producer/allocation changes, shader compilation, dispatch, frame or full qualification.'}, indent=2) + '\n')
print(json.dumps({family: {key: value for key, value in result.items() if key != 'matches'}
                  for family, result in results.items()}, indent=2))
