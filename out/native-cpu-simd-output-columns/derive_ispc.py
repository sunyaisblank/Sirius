"""Translate the qualified scalar headers to per-row SPMD source, not physics.

Default checks existing derived bytes. --write creates missing outputs only.
The original 9380 literal instruction calls and every owned-row body are kept.
"""
from pathlib import Path
import hashlib
import json
import re
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
basis = root / 'bin/linux-gcc/native-cpu-static-program/typed-batch'
destination = root / 'bin/linux-gcc/native-cpu-simd-output-columns/ispc-source'
assert sys.argv[1:] in ([], ['--write'])
write = sys.argv[1:] == ['--write']


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert write, 'Missing derived source: ' + str(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


def end_of_block(text, opening):
    depth, end = 1, opening + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return end


def close_delimiter(text, opening, left, right):
    depth, end = 1, opening + 1
    while depth:
        depth += (text[end] == left) - (text[end] == right)
        end += 1
    return end - 1


def uniform_register_access(text):
    """Clone translated private-register helpers; preserve generic bodies."""
    before = text
    names = {name + '_0': name + 'Uniform_0' for name in
             ('ScratchAddress', 'ScratchRead', 'ScratchWrite', 'ReadValue', 'WriteValue')}
    insertions = []
    order = []
    for old, new in names.items():
        match, = re.finditer(r'^static (?:uint32_t|RetainedTriple_0|void) ' + old +
                             r'\([^\n]+\)\n\{', text, flags=re.M)
        opening = match.end() - 1
        end = end_of_block(text, opening)
        original = text[match.start():end]
        header, body = original.split('{', 1)
        header, qualified = re.subn(r'\buint32_t ((?:index|component|prefix)_\d+)\b',
                                    r'uniform uint32_t \1', header)
        assert qualified == (3 if old.startswith('Scratch') else 1)
        if old == 'ScratchAddress_0':
            assert not re.search(r'\browBase_\d+\b', body)
            header = header.replace('static uint32_t ', 'static uniform uint32_t ', 1)
        header = header.replace(old + '(', new + '(', 1)
        clone = header + '{' + body
        for source, target in names.items():
            clone = re.sub(r'\b' + source + r'(?=\()', target, clone)
        # Reverse only approved names/qualifiers to prove exact body cloning.
        restored = clone
        for source, target in names.items():
            restored = re.sub(r'\b' + target + r'(?=\()', source, restored)
        restored_header, restored_body = restored.split('{', 1)
        restored_header = restored_header.replace('uniform uint32_t ', 'uint32_t ')
        assert restored_header + '{' + restored_body == original
        addition = '\n\n// Uniform canonical register access; original generic helper retained.\n' + clone
        insertions.append((end, addition))
        order.append((match.start(), new))
    for end, addition in sorted(insertions, reverse=True):
        text = text[:end] + addition + text[end:]

    match, = re.finditer(r'^static bool EvaluateOneStatic_0\(([^\n]+)\)\n\{', text, flags=re.M)
    opening = match.end() - 1
    end = end_of_block(text, opening)
    original_static = text[match.start():end]
    header, body = original_static.split('{', 1)
    header, qualified = re.subn(r'\buint32_t ((?:destination|a|b|c)_\d+)\b',
                                r'uniform uint32_t \1', header)
    assert qualified == 4
    for name, expected in (('ReadValue_0', 3), ('WriteValue_0', 1)):
        body, changed = re.subn(r'\b' + name + r'(?=\()', names[name], body)
        assert changed == expected
    revised_static = header + '{' + body
    text = text[:match.start()] + revised_static + text[end:]
    restored = text.replace(revised_static, original_static, 1)
    for _, addition in insertions:
        assert restored.count(addition) == 1
        restored = restored.replace(addition, '', 1)
    assert restored == before
    return text, [name for _, name in sorted(order)]


def relative_output_columns(text, slug, row_words):
    """Change only private output location and the equivalent full clear."""
    before = text
    field = '    uint32_t simd_output_row_base;\n'
    assert text.count(field) == 1
    text = text.replace(field, '', 1)
    declaration, = re.finditer(r'    uint32_t (rowBase_\d+) = (row_\d+) \* (rowWords_\d+);', text)
    original_base = declaration[0]
    base, real_row, words = declaration.groups()
    new_base = f'    uniform uint32_t {base} = 0U;'
    text = text.replace(original_base, new_base, 1)
    assert f'uint32_t {words} = {row_words}U;' in text
    parameter_count = 0
    signatures = list(re.finditer(r'^static [^\n]+\([^\n]+\)\n\{', text, flags=re.M))
    for match in reversed(signatures):
        header = match[0]
        changed, count = re.subn(r'\buint32_t (rowBase_\d+)\b', r'uniform uint32_t \1', header)
        assert count <= 1
        parameter_count += count
        text = text[:match.start()] + changed + text[match.end():]
    assert parameter_count == {'transport': 173, 'endpoint': 212, 'dopri_phase': 58}[slug]

    main_start = text.index('static void ComputeOwnedRow_0(')
    main_open = text.index('{', main_start)
    main_end = end_of_block(text, main_open)
    clear, = re.finditer(r'\*\(&OUTPUT_WORD\((\w+), ' + base + r' \+ (i_\d+)\)\) = 0U;', text[main_start:main_end])
    context, index = clear.groups()
    clear_position = main_start + clear.start()
    loop_start = text.rfind('    for(;;)\n    {', main_start, clear_position)
    assert loop_start >= main_start
    loop_end = end_of_block(text, text.index('{', loop_start))
    original_clear = text[loop_start:loop_end]
    normalized = re.sub(r'\s+', '', re.sub(r'^#line[^\n]*\n', '', original_clear, flags=re.M))
    expected = (f'for(;;){{if({index}<{words}){{}}else{{break;}}'
                f'*(&OUTPUT_WORD({context},{base}+{index}))=0U;{index}={index}+1U;}}')
    assert normalized == expected
    assert re.search(r'uint32_t ' + index + r' = lane_\d+;', text[main_start:loop_start])
    after = text[loop_end:main_end]
    first_use = re.search(r'\b' + index + r'\b', after)
    reset = re.search(r'^    ' + index + r' = 0U;', after, flags=re.M)
    assert first_use and reset and first_use.start() == reset.start() + 4
    # The old clear counter is never read before this unconditional reset.
    new_clear = (f'    for (uniform uint32_t clear_column = 0U; clear_column < {row_words}U; ++clear_column) {{\n'
                 f'        *(&OUTPUT_WORD({context}, {base} + clear_column)) = 0U;\n'
                 '    }')
    text = text[:loop_start] + new_clear + text[loop_end:]
    restored = text.replace(new_clear, original_clear, 1)
    restored = re.sub(r'\buniform uint32_t (rowBase_\d+)\b', r'uint32_t \1', restored)
    restored = restored.replace(new_base.replace('uniform ', ''), original_base, 1)
    restored = restored.replace('struct KernelContext_0\n{\n', 'struct KernelContext_0\n{\n' + field, 1)
    assert restored == before, 'Only approved output-address/clear changes are permitted'
    return text, {'uniform_forwarded_row_base_parameters': parameter_count,
                  'top_row_base_zeroed': 1, 'context_absolute_base_fields_removed': 1,
                  'complete_clear_loops_replaced': 1,
                  'old_clear_counter_first_later_use': reset[0].strip(),
                  'real_input_row_unchanged': real_row,
                  'clear_words_in_original_order': row_words,
                  'all_other_body_bytes_unchanged_after_inverse': True}


def remove_redundant_copy_in(entry, row_words):
    before = entry
    start = entry.index('    // Preserve even original unmodified words')
    end = entry.index('    for (uniform uint32_t invocation', start)
    copy_in = entry[start:end]
    expected = f'''    // Preserve even original unmodified words on every active row. Only the
    // final copy publishes the tile; masked-off rows are never accessed.
    if (programIndex < groups) {{
        for (uniform uint32_t column = 0; column < {row_words}U; ++column)
            staged_output[column * programCount + programIndex] =
                output[(first_x + programIndex) * {row_words}U + column];
    }}
'''
    assert copy_in == expected
    replacement = '    // Full admission proves every active owned row reaches its complete clear.\n'
    entry = entry[:start] + replacement + entry[end:]
    assignment = f'            context.simd_output_row_base = (first_x + programIndex) * {row_words}U;\n'
    assert entry.count(assignment) == 1
    entry = entry.replace(assignment, '', 1)
    restored = entry.replace(replacement, copy_in, 1)
    restored = restored.replace('            context.globalParams_0 = &params;\n',
                                '            context.globalParams_0 = &params;\n' + assignment, 1)
    assert restored == before
    return entry


PRELUDE = r'''
typedef uint32 uint32_t;
typedef int32 int32_t;
struct RowID { uint32_t x; uint32_t y; uint32_t z; };
struct WordDimensions { uint32_t x; uint32_t y; };
static WordDimensions uint2(uint32_t x, uint32_t y) {
    WordDimensions r; r.x=x; r.y=y; return r;
}
struct WordBuffer {
    uniform uint32_t * uniform data;
    uniform uint64 count;
};
struct GlobalParams_0 {
    uniform WordBuffer inputs_0;
    uniform WordBuffer outputs_0;
};
// Original Slang/GCC prelude uses libm fmaxf. Preserve NaN and zero-sign
// behavior through the strict C++ libm adapter for exceptional/tied zeros.
// A unique symbol prevents ISPC from replacing this call as a math builtin.
extern "C" uniform float sirius_simd_fmaxf(uniform float a, uniform float b);
static float F32_max(float a, float b) {
    float r = a > b ? a : b;
    bool edge = isnan(a) || isnan(b) ||
        (((intbits(a) & 0x7fffffffU) == 0U) &&
         ((intbits(a) ^ intbits(b)) == 0x80000000U));
    if (edge) {
        foreach_active (i) {
            r = sirius_simd_fmaxf(extract(a, i), extract(b, i));
        }
    }
    return r;
}
static uint32_t F32_asuint(float x) { return intbits(x); }
static float U32_asfloat(uint32_t x) { return floatbits(x); }
static float F32_sqrt(float x) { return sqrt(x); }
static uint32_t U32_firstbithigh(uint32_t x) {
    if (x == 0U) return 0xffffffffU;
    return 31U - count_leading_zeros(x);
}
#define LOAD_WORD(buffer, index) ((buffer).data[(index)])
#define OUTPUT_WORD(context, index) ((context)->globalParams_0->outputs_0.data[ \
    ((index) - (context)->simd_output_row_base) * programCount + programIndex])
'''

G_PRELUDE = PRELUDE
old_output_offset = '((index) - (context)->simd_output_row_base)'
assert PRELUDE.count(old_output_offset) == 1
PRELUDE = PRELUDE.replace(old_output_offset, '(index)', 1)
g_generator = root / 'out/native-cpu-simd-rows/derive_ispc.py'
assert pin(g_generator.read_bytes())['sha256'] == '8875e5bb3978e7c3072a3d053cc3ab215368d07995411be22cb4676573363722'
g_derivation = json.loads((g_generator.parent / 'ispc-uniform-derivation.json').read_text())

authority = json.loads((root / 'out/native-cpu-static-program/header-preparation.json').read_text())
records = {}
total = 0
for slug in ('transport', 'endpoint', 'dopri_phase'):
    source = basis / (slug + '.native.hpp')
    raw = source.read_text()
    assert pin(source.read_bytes()) == authority['headers'][slug]['derived']
    # The original dispatcher group has 64 invocations, with an owned body in
    # lane zero. The export below keeps that logical loop for each active row.
    start = raw.index('namespace sirius_generated_' + slug + ' {')
    start = raw.index('\n', start) + 1
    end = raw.index('void sirius_' + slug + '_ComputeMain(void*')
    text = raw[start:end]
    text = re.sub(r'#ifdef SLANG_PRELUDE_NAMESPACE\nusing namespace SLANG_PRELUDE_NAMESPACE;\n#endif\n', '', text)
    # GlobalParams layout is local to this translation, never shared as a C ABI
    # object. Its storage pointers/counts are uniform across independent rows.
    gp = text.index('struct GlobalParams_0\n{')
    gp_end = end_of_block(text, text.index('{', gp))
    assert text[gp_end:gp_end + 1] == ';'
    text = text[:gp] + text[gp_end + 1:]
    text = text.replace('struct KernelContext_0\n{',
                        'struct KernelContext_0\n{\n    uint32_t simd_output_row_base;')
    text = text.replace('GlobalParams_0* globalParams_0;', 'uniform GlobalParams_0* uniform globalParams_0;')
    counts = {}
    text, counts['private_pointers'] = re.subn(
        r'\b(RetainedTriple_0|KernelContext_0|float|uint32_t)\s*\*',
        r'varying \1 * uniform', text)
    # Two-dimensional immutable Butcher arrays copy their original literal rows.
    rows = {}
    def flat_array(match):
        kind, number, name, initializer = match.groups()
        rows[name] = initializer
        return f'static uniform const {kind} {name}[{number}] = {initializer};'
    text, counts['constant_arrays'] = re.subn(
        r'static const FixedArray<(int32_t), (\d+)>\s+(\w+) = (\{[^;]+\});', flat_array, text)
    def matrix(match):
        kind, columns, nrows, name, contents = match.groups()
        names = [n.strip() for n in contents.split(',')]
        assert len(names) == int(nrows) and all(n in rows for n in names)
        initializer = '{' + ','.join(rows[n] for n in names) + '}'
        return f'static uniform const {kind} {name}[{nrows}][{columns}] = {initializer};'
    text, counts['constant_matrices'] = re.subn(
        r'static const FixedArray<FixedArray<(int32_t), (\d+)> , (\d+)>\s+(\w+) = \{([^;]+)\};', matrix, text)
    text, counts['private_arrays'] = re.subn(
        r'FixedArray<(uint32_t), (\d+)>\s+(\w+);', r'\1 \3[\2];', text)
    text = re.sub(r'Vector<uint32_t, 3>\s*', 'RowID ', text)
    text = re.sub(r'Vector<uint32_t, 2>\s*', 'WordDimensions ', text)
    # Methods perform the same uint word load. Their original indices and
    # runtime input reads remain untouched, including the E projection handoff.
    member = r'\w+(?:->\w+)*->(?:inputs_0|outputs_0)'
    # Match either a complete parenthesized receiver or a bare receiver.
    # An optional opening parenthesis would also consume a caller's '(' .
    buffer = r'(\(' + member + r'\)|' + member + r')'
    text, counts['dimension_reads'] = re.subn(
        buffer + r'\.GetDimensions\(&([\w]+), &([\w]+)\);',
        r'\2 = (\1).count; \3 = 4U;', text)
    assert counts['dimension_reads'] == 2
    # Public output words use a private SoA tile during execution. Each source
    # index remains an absolute original output index; subtraction only maps it
    # to that row's column. This changes storage, never an arithmetic expression.
    for method, left, right in (('word_loads', '(', ')'), ('word_indices', '[', ']')):
        pattern = re.compile(buffer + (r'\.Load\(' if method == 'word_loads' else r'\s*\['))
        offset = 0
        count = 0
        while match := pattern.search(text, offset):
            opening = match.end() - 1
            closing = close_delimiter(text, opening, left, right)
            receiver = match[1].strip('()')
            index = text[opening + 1:closing]
            context = receiver.split('->')[0]
            if receiver.endswith('outputs_0'):
                replacement = f'OUTPUT_WORD({context}, {index})'
            elif method == 'word_loads':
                replacement = f'LOAD_WORD({receiver}, {index})'
            else:
                replacement = f'({receiver}).data[{index}]'
            text = text[:match.start()] + replacement + text[closing + 1:]
            offset = match.start() + len(replacement)
            count += 1
        counts[method] = count
    row_words = {'transport': 4699, 'endpoint': 3829, 'dopri_phase': 4149}[slug]
    # Admission has already established exactly this output stride. Keep the
    # original validation; specialize only its established quotient, allowing
    # absolute-to-private output offsets to cancel without integer divisions.
    text, counts['admitted_output_stride'] = re.subn(
        r'(uint32_t rowWords_\d+ = )no_\d+ / rows_\d+;',
        r'\g<1>' + str(row_words) + 'U;', text)
    assert counts['admitted_output_stride'] == 1
    text, counts['compile_assertions'] = re.subn(r'^\s*static_assert\([^;]+;\n', '\n', text, flags=re.M)
    text = text.replace('constexpr uint32_t op_0 = StaticOp;', 'uniform const uint32_t op_0 = StaticOp;')
    # ISPC supports the same function-template/integer-control constructs, but
    # scalar conversion uses C cast spelling instead of C++ constructor syntax.
    for kind in ('int32_t', 'uint32_t', 'float', 'int'):
        text = re.sub(r'\b' + kind + r'\(', '(' + kind + ')(', text)
    text, uniform_helper_order = uniform_register_access(text)
    counts['uniform_register_helper_clones'] = len(uniform_helper_order)
    counts['uniform_canonical_operand_parameters'] = 4
    g_text = text
    text, output_columns = relative_output_columns(text, slug, row_words)
    counts['uniform_forwarded_output_bases'] = output_columns['uniform_forwarded_row_base_parameters']
    counts['zero_output_base_locals'] = 1
    counts['removed_context_output_base_fields'] = 1
    counts['equivalent_uniform_full_output_clear_loops'] = 1
    counts['removed_export_output_base_assignments'] = 1
    counts['removed_initial_global_output_copy_loops'] = 1
    calls = re.findall(r'EvaluateOneStatic_0<(\d+)U>\(([^\n]+)\)', text)
    original_calls = re.findall(r'EvaluateOneStatic_0<(\d+)U>\(([^\n]+)\)', raw)
    assert calls == original_calls
    assert re.findall(r'^.*EvaluateOneStatic_0<\d+U>.*$', text, re.M) == re.findall(r'^.*EvaluateOneStatic_0<\d+U>.*$', raw, re.M)
    expected = sum(p['instructions'] for p in authority['headers'][slug]['programs'])
    assert len(calls) == expected
    total += len(calls)
    assert not any(s in text for s in ('FixedArray<', 'Vector<', '.Load(', '.GetDimensions(', 'slang_bit_cast', 'static_assert', 'constexpr'))
    entry = f'''
export void sirius_simd_{slug}_packet(
    uniform uint32_t * uniform input, uniform uint32_t input_words,
    uniform uint32_t * uniform output, uniform uint32_t output_words,
    uniform uint32_t first_x, uniform uint32_t groups) {{
    uniform GlobalParams_0 params;
    params.inputs_0.data = input; params.inputs_0.count = input_words;
    uniform uint32_t staged_output[{row_words} * programCount];
    params.outputs_0.data = staged_output; params.outputs_0.count = output_words;
    // Preserve even original unmodified words on every active row. Only the
    // final copy publishes the tile; masked-off rows are never accessed.
    if (programIndex < groups) {{
        for (uniform uint32_t column = 0; column < {row_words}U; ++column)
            staged_output[column * programCount + programIndex] =
                output[(first_x + programIndex) * {row_words}U + column];
    }}
    for (uniform uint32_t invocation = 0; invocation < 64U; ++invocation) {{
        if (invocation == 0U) {{
            varying KernelContext_0 context;
            memset((void * uniform)&context, 0, sizeof(context));
            context.globalParams_0 = &params;
            context.simd_output_row_base = (first_x + programIndex) * {row_words}U;
            RowID id; id.x = first_x + programIndex; id.y = 0U; id.z = 0U;
            if (programIndex < groups) ComputeOwnedRow_0(id, invocation, &context);
        }}
    }}
    if (programIndex < groups) {{
        for (uniform uint32_t column = 0; column < {row_words}U; ++column)
            output[(first_x + programIndex) * {row_words}U + column] =
                staged_output[column * programCount + programIndex];
    }}
}}
'''
    assert pin((G_PRELUDE + g_text + entry).encode()) == g_derivation['stages'][slug]['derived']
    entry = remove_redundant_copy_in(entry, row_words)
    derived = PRELUDE + text + entry
    assert 'simd_output_row_base' not in derived
    output = destination / (slug + '.ispc')
    publish(output, derived.encode())
    records[slug] = {'basis_path': str(source.relative_to(root)), 'basis': pin(source.read_bytes()),
        'derived_path': str(output.relative_to(root)), 'derived': pin(derived.encode()),
        'approved_syntax_replacement_counts': counts, 'literal_instruction_calls_preserved': len(calls),
        'uniform_register_helper_order': uniform_helper_order,
        'relative_output_columns': output_columns,
        'private_output_layout': {'row_words': row_words,
            'index_mapping': 'row_relative_column * programCount + programIndex',
            'active_rows_copied_in': False, 'active_rows_copied_out': True,
            'complete_original_owned_clear_retained': True, 'inactive_rows_accessed': False,
            'all_original_output_values_preserved': True},
        'preclear_admission_domain': {
            'authority': 'out/native-cpu-active-prefix/native_packet.hpp',
            'maximum_capacity': 65536,
            'minimum_program_header_words': {'transport': 43, 'endpoint': 103, 'dopri_phase': 203}[slug],
            'actual_canonical_program_words': {'transport': 15967, 'endpoint': 21143, 'dopri_phase': 10689}[slug],
            'minimum_output_prefix_words': {'transport': 2404, 'endpoint': 769, 'dopri_phase': 1004}[slug],
            'conditions': 'Full B admission and immutable views prove positive capacity<=65536, original input/program extents, output capacity*stride/divisibility, active rows inside capacity and sufficient output prefix. Export supplies y=z=0 and owned invocation0. Thus every active row fully clears before any numerical/input-field rejection; no previous global output word is read or needed. Scalar fallback remains unchanged.'},
        'logical_invocations_per_active_row': 64, 'physical_math_rederived': False,
        'limits': 'Source derivation only. Uniform/global versus varying/private storage adaptation and prelude intrinsics need independent review and actual compiled/numerical proof. No object/execution/benefit/frame claim.'}
assert total == 9380
publish(work / 'ispc-derivation.json', (json.dumps({'instructions': total, 'stages': records,
    'compiler_user_guide': 'https://ispc.github.io/ispc.html',
    'scope': 'Complete T/E/DP plus all5canonicalprograms. Original arithmetic/control expressions copied with enumerated C-like syntax/type/storage adaptations; original literal call order retained. No field/trial/RHS/quartic/reference/precision/workload/count/budget change. Not compiled or executed.'}, indent=2) + '\n').encode())
print(json.dumps({'derived_stages': 3, 'instructions': total, 'compiled_or_executed': False}))
