"""Derive scalar row-control bridges and field-SIMD layer executors for all3 stages.

Existing scalar arithmetic, RK/control/refusal/publication and all canonical
instructions remain. Only the known-layer dispatch delegates to the new ISPC
executor. No compile, link, provider mutation or physics occurs here.
"""
from pathlib import Path
import hashlib
import json
import re
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-field-simd'
assert sys.argv[1:] in ([], ['--write'])


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write'], path
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)


def function(text, name):
    pattern = r'(?m)^static (?:inline )?[^\n{;]+\b' + re.escape(name) + r'\('
    match = re.search(pattern, text)
    assert match, name
    opening = text.index('{', match.end()); end = opening + 1; depth = 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}'); end += 1
    return match.start(), opening, end


schedule = json.loads((W / 'schedule-preparation.json').read_text())
base = json.loads((R / 'out/native-cpu-static-program/header-preparation.json').read_text())
assert schedule['instructions'] == base['instructions'] == 9380
record = {'scope': 'Source derivation only. Scalar original row control calls field-SIMD only at complete canonical layer boundaries; original arithmetic, coefficients, roots, publication, generic fallback and logical64 calls retained. No runtime or benefit/scientific qualification.',
          'source_inputs': {}, 'stages': {}, 'new_CPP_compilations_planned': 3,
          'new_ISPC_compilations_planned': 3, 'field_width': 8}
for source in (W / 'prepare_schedules.py', W / 'schedule-preparation.json', Path(__file__).resolve(),
               R / 'out/native-cpu-static-program/header-preparation.json'):
    record['source_inputs'][str(source.relative_to(R))] = pin(source.read_bytes())

for slug in ('transport', 'endpoint', 'dopri_phase'):
    scalar_path = R / base['headers'][slug]['derived_path']
    vector_path = R / ('bin/linux-gcc/native-cpu-zero-term/ispc-source/' + slug + '.ispc')
    wrapper_path = R / ('bin/linux-gcc/native-cpu-slang-bridge/typed-batch/' + slug + '.cpp')
    table_path = R / base['headers'][slug]['table_path']
    packed_path = R / schedule['families'][slug]['generated_path']
    scalar, vector, wrapper = scalar_path.read_text(), vector_path.read_text(), wrapper_path.read_text()
    assert pin(scalar.encode()) == base['headers'][slug]['derived']
    assert pin(packed_path.read_bytes()) == schedule['families'][slug]['generated']
    scratch = int(re.search(r'retainedScratch_0\[(\d+)\]', vector)[1])
    coefficient = re.search(r'retainedCoefficients_0\[(\d+)\]', vector)
    coefficients = int(coefficient[1]) if coefficient else 0
    export = 'sirius_simd_' + slug + '_field_layer'
    cpp_namespace = 'sirius_generated_' + slug
    new_namespace = cpp_namespace + '_field'
    known_start, known_open, known_end = function(scalar, 'EvaluateKnownLayer_0')
    known_original = scalar[known_start:known_end]
    assert re.search(r'uint32_t programBase, uint32_t count, uint32_t registers,', known_original)
    declaration = ('extern "C" bool ' + export + '(const uint32_t*, uint32_t, uint32_t*, uint32_t,\n'
        '    uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t,\n'
        '    uint32_t*, uint32_t*, uint32_t*) noexcept;\n\n')
    coefficient_pointer = '&context->retainedCoefficients_0[0]' if coefficients else 'nullptr'
    bridge = '{\n    return ' + export + '(\n'
    bridge += '        context->globalParams_0->inputs_0.data, static_cast<uint32_t>(context->globalParams_0->inputs_0.count),\n'
    bridge += '        context->globalParams_0->outputs_0.data, static_cast<uint32_t>(context->globalParams_0->outputs_0.count),\n'
    bridge += '        programBase, count, registers, layer, start, end, row, rowBase,\n'
    bridge += '        &context->retainedScratch_0[0], ' + coefficient_pointer + ', &context->retainedStatus_0);\n}'
    derived_scalar = scalar[:known_open] + bridge + scalar[known_end:]
    namespace_marker = 'namespace ' + cpp_namespace + ' {'
    assert derived_scalar.count(namespace_marker) == 1
    derived_scalar = derived_scalar.replace(namespace_marker, declaration + 'namespace ' + new_namespace + ' {', 1)
    # Slang exports use C linkage; namespace renaming alone cannot separate them.
    # Only the generated C thread/group wrappers need a new C-linkage label.
    # The similarly named internal entry already gains the new C++ namespace.
    labels = re.findall(r'void (sirius_\w+ComputeMain)_Thread\(', scalar)
    labels = sorted(set(labels)); assert len(labels) == 1, labels
    entry_prefix = labels[0]
    new_entry_prefix = entry_prefix.replace('ComputeMain', 'FieldComputeMain')
    derived_scalar = derived_scalar.replace(entry_prefix, new_entry_prefix)
    inverse = derived_scalar.replace(new_entry_prefix, entry_prefix).replace(
        declaration + 'namespace ' + new_namespace + ' {', namespace_marker, 1)
    inverse_start, inverse_open, inverse_end = function(inverse, 'EvaluateKnownLayer_0')
    inverse = inverse[:inverse_open] + scalar[known_open:known_end] + inverse[inverse_end:]
    assert inverse == scalar
    derived_wrapper = wrapper.replace(cpp_namespace, new_namespace).replace(entry_prefix, new_entry_prefix)
    original_packet = 'sirius_native_' + slug + '_try_packet'
    field_packet = 'sirius_native_' + slug + '_field_try_packet'
    assert derived_wrapper.count(original_packet + '(') == 1
    derived_wrapper = derived_wrapper.replace(original_packet + '(', field_packet + '(', 1)
    derived_wrapper += '\nstatic_assert(std::is_same_v<decltype(&' + field_packet + '), decltype(&' + original_packet + ')>);\n'
    cpp_header = B / 'cpp' / (slug + '.native.hpp')
    cpp_source = B / 'cpp' / (slug + '.cpp')
    publish(cpp_header, derived_scalar.encode()); publish(cpp_source, derived_wrapper.encode())

    # Keep the exact S arithmetic functions. The shared scratch bank is the
    # scalar caller's existing array; distinct packed destinations are proven.
    matches = list(re.finditer(r'(?m)^static (?:inline )?[^\n{;]+\b(\w+)\(', vector))
    functions = {}
    for match in matches:
        name = match[1]
        start, opening, end = function(vector, name)
        assert name not in functions
        functions[name] = (start, end, vector[start:end])
    start, _opening, end = function(vector, 'EvaluateOneStatic_0')
    field_one = vector[start:end].replace('EvaluateOneStatic_0(', 'EvaluateOneField_0(', 1)
    field_one = field_one.replace('ReadValueUniform_0(', 'ReadValue_0(').replace('WriteValueUniform_0(', 'WriteValue_0(')
    signature_end = field_one.index('{')
    signature = field_one[:signature_end]
    parameters = re.findall(r'uniform uint32_t (\w+)', signature)
    assert len(parameters) >= 5
    for name in parameters[:4]:
        signature = signature.replace('uniform uint32_t ' + name, 'uint32_t ' + name, 1)
    field_one = 'template<uint32_t StaticOp>\n' + signature + field_one[signature_end:]
    pending = set(re.findall(r'\b(\w+)\s*\(', field_one)) & functions.keys()
    closure = set()
    while pending:
        name = pending.pop()
        if name in closure:
            continue
        closure.add(name)
        pending.update((set(re.findall(r'\b(\w+)\s*\(', functions[name][2])) & functions.keys()) - closure)
    assert 'EvaluateOneStatic_0' not in closure and 'ComputeOwnedRow_0' not in closure
    assert not any('CacheCoefficient' in name for name in closure)
    first_line = vector.index('\n#line')
    type_start = vector.index('struct RetainedTriple_0')
    type_end = vector.index('};', type_start) + 2
    prelude = vector[:type_start] + vector[type_start:type_end]
    macro = '#define OUTPUT_WORD(context, index) ((context)->globalParams_0->outputs_0.data[ \\\n    (index) * programCount + programIndex])'
    assert prelude.count(macro) == 1
    prelude = prelude.replace(macro, '#define OUTPUT_WORD(context, index) ((context)->globalParams_0->outputs_0.data[(index)])')
    context = '''
struct KernelContext_0 {
    uniform GlobalParams_0 * uniform globalParams_0;
    uniform uint32_t * uniform retainedScratch_0;
    uint32_t retainedStatus_0;
    uniform uint32_t * uniform retainedCoefficients_0;
};
'''
    helpers = '\n'.join(functions[name][2] for name in sorted(closure, key=lambda name: functions[name][0]) if functions[name][0] >= first_line)
    field_source = prelude + '\n' + context + '\n' + helpers + '\n' + field_one + '\n'
    field_source += '#include "' + str(packed_path) + '"\n'
    words = [int(value) for value in re.findall(r'(\d+)u', table_path.read_text())]
    field_source += 'static const uniform uint32_t SiriusFieldWire[' + str(len(words)) + '] = {\n'
    field_source += ',\n'.join('    ' + ', '.join(str(value) + 'U' for value in words[i:i+16]) for i in range(0,len(words),16)) + '\n};\n'
    field_source += '\nexport uniform bool ' + export + '''(
    uniform uint32_t * uniform input, uniform uint32_t input_words,
    uniform uint32_t * uniform output, uniform uint32_t output_words,
    uniform uint32_t program_base, uniform uint32_t count, uniform uint32_t registers,
    uniform uint32_t layer, uniform uint32_t start, uniform uint32_t end,
    uniform uint32_t row, uniform uint32_t row_base,
    uniform uint32_t * uniform scratch, uniform uint32_t * uniform coefficients,
    uniform uint32_t * uniform status) {
    if (input == NULL || output == NULL || scratch == NULL || status == NULL || input_words == 0U)
        return false;
'''
    if coefficients:
        field_source += '    if (coefficients == NULL) return false;\n'
    stride = base['headers'][slug]['input_row_words_from_original_packet_guard']
    field_source += f'    uniform uint32_t table_base = 1U + input[0] * {stride}U;\n'
    for family, item in enumerate(schedule['families'][slug]['programs']):
        prefix = f'SiriusField_{slug}_{family}_8'
        instruction_base = item['base'] + 3 + item['roots']
        field_source += f'''    if (program_base == table_base + {item['base']}U && count == {item['instructions']}U && registers == {item['registers']}U && layer < {item['layers']}U) {{
        uniform uint32_t offsets_base = {instruction_base + 5*item['instructions']}U;
        if (start != SiriusFieldWire[offsets_base + layer] || end != SiriusFieldWire[offsets_base + layer + 1U]) return false;
        uniform GlobalParams_0 params;
        params.inputs_0.data = input; params.inputs_0.count = input_words;
        params.outputs_0.data = output; params.outputs_0.count = output_words;
        varying KernelContext_0 context;
        context.globalParams_0 = &params;
        context.retainedScratch_0 = scratch;
        context.retainedCoefficients_0 = coefficients;
        context.retainedStatus_0 = *status;
        for (uniform uint32_t pack = {prefix}_LayerPacks[layer]; pack < {prefix}_LayerPacks[layer+1U]; ++pack) {{
            uniform uint32_t op = {prefix}_Opcode[pack];
            if (programIndex < {prefix}_ActiveNodes[pack]) {{
                uint32_t node = {prefix}_Instruction[pack*8U + programIndex];
                uint32_t address = {instruction_base}U + node*5U;
                uint32_t dst = SiriusFieldWire[address+1U];
                uint32_t a = SiriusFieldWire[address+2U];
                uint32_t b = SiriusFieldWire[address+3U];
                uint32_t c = SiriusFieldWire[address+4U];
                bool handled = false;
                switch (op) {{
'''
        for op in range(12):
            field_source += f'                case {op}U: handled = EvaluateOneField_0<{op}U>(dst,a,b,c,row,row_base,registers,&context); break;\n'
        field_source += '''                default: break;
                }
                if (!handled) context.retainedStatus_0 = 1U;
            }
        }
        if (any(context.retainedStatus_0 != 0U)) *status = 1U;
        return true;
    }
'''
    field_source += '    return false;\n}\n'
    ispc_source = B / 'ispc' / (slug + '.ispc')
    publish(ispc_source, field_source.encode())
    for path in (scalar_path, vector_path, wrapper_path, table_path, packed_path):
        record['source_inputs'][str(path.relative_to(R))] = pin(path.read_bytes())
    record['stages'][slug] = {'cpp_header': {'path': str(cpp_header.relative_to(R)), **pin(derived_scalar.encode())},
        'cpp_source': {'path': str(cpp_source.relative_to(R)), **pin(derived_wrapper.encode())},
        'ispc_source': {'path': str(ispc_source.relative_to(R)), **pin(field_source.encode())},
        'original_known_layer': pin(known_original.encode()), 'scalar_inverse_restores_original_header': True,
        'arithmetic_helpers_copied_verbatim': {name: pin(functions[name][2].encode()) for name in sorted(closure)},
        'scratch_words': scratch, 'coefficient_words': coefficients, 'export': export,
        'field_packet_export': field_packet, 'original_scalar_fallback_preserved_separately': True,
        'program_instruction_count': sum(item['instructions'] for item in schedule['families'][slug]['programs'])}
assert sum(item['program_instruction_count'] for item in record['stages'].values()) == 9380
publish(W / 'executor-derivation.json', (json.dumps(record, indent=2) + '\n').encode())
print(json.dumps({'derived_stages': 3, 'CPP_control_bridges': 3, 'ISPC_layer_executors': 3,
                  'canonical_instructions': 9380, 'compiler_or_workload_executed': False}))
