"""Derive only pinned private canonical metadata; never compile or execute it."""
from pathlib import Path
import hashlib
import json
import re
import sys

R = Path(__file__).resolve().parents[2]
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-canonical-metadata'
STAGES = ('dopri_phase', 'transport', 'endpoint')
PROFILES = {
    'dopri_phase': dict(stride=1810, roots=200, row='row_2', rowbase='rowBase_6',
        base='programBase_1', count='count_0', registers='registers_1', lane='lane_0',
        context='kernelContext_8', layers=None, programs=[(0, 2088, 629, 45, None)]),
    'transport': dict(stride=230, roots=40, row='row_3', rowbase='rowBase_7',
        base='programBase_1', count='count_0', registers='registers_1', lane='lane_1',
        context='kernelContext_12', layers='layers_0', programs=[(0, 2587, 459, 117, None), (13096, 557, 459, 42, None)]),
    'endpoint': dict(stride=225, roots=100, row='row_1', rowbase='rowBase_10',
        base='programBase_2', count='count_1', registers='registers_2', lane='lane_0',
        context='kernelContext_12', layers='layers_1', programs=[(0, 2205, 612, 111, 250), (11240, 1943, 612, 84, 202)]),
}
BEGIN = '    /* BEGIN pinned canonical metadata specialization. */\n'
END = '    /* END pinned canonical metadata specialization; exact fallback follows. */\n'


def pin_bytes(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    return pin_bytes(path.read_bytes())


def function(source, name):
    match = re.search(r'^static [^\n]*\b' + re.escape(name) + r'\([^;{]*\)\n\{', source, re.M)
    assert match, name
    opening = match.end() - 1
    depth = 1
    i = opening + 1
    while depth:
        assert i < len(source), name
        depth += (source[i] == '{') - (source[i] == '}')
        i += 1
    return source[match.start():i], opening + 1


def program_words(stage):
    text = (R / f'bin/linux-gcc/native-cpu-slang-bridge/typed-batch/{stage}.program.hpp').read_text()
    match = re.search(r'inline constexpr std::uint32_t ' + stage + r'_program\[\] = \{(.*?)\};', text, re.S)
    assert match
    assert not re.sub(r'\d+u|[\s,]', '', match[1])
    return [int(n) for n in re.findall(r'(\d+)u', match[1])]


def canonical_programs(source, stage):
    p = PROFILES[stage]
    words = program_words(stage)
    known, _ = function(source, 'EvaluateKnownLayer_0')
    groups = re.findall(r'if\(programBase == tableBase \+ (\d+)U && count == (\d+)U && registers == (\d+)U\)\n'
        r'    \{\n        switch\(layer\)\n        \{(.*?)        default: return false;\n        \}\n    \}', known, re.S)
    assert len(groups) == len(p['programs'])
    result = []
    for expected, (offset, count, registers, cases) in zip(p['programs'], groups):
        offset, count, registers = map(int, (offset, count, registers))
        assert (offset, count, registers) == expected[:3]
        layers = re.findall(r'case (\d+)U:\n            if\(start != (\d+)U \|\| end != (\d+)U\) return false;\n'
            r'            (EvaluateStaticLayer_\d+_\d+_0)\(row, rowBase, context\);\n            return true;', cases)
        assert len(layers) == expected[3]
        assert words[offset:offset + 3] == [count, registers, p['roots']]
        roots = words[offset + 3:offset + 3 + p['roots']]
        assert len(roots) == p['roots'] and all(index < registers for index in roots)
        instruction_base = offset + 3 + p['roots']
        layer_base = instruction_base + 5 * count
        boundaries = words[layer_base:layer_base + len(layers) + 1]
        assert len(boundaries) == len(layers) + 1 and boundaries[0] == 0 and boundaries[-1] == count
        literal = []
        for number, (layer, start, end, name) in enumerate(layers):
            layer, start, end = map(int, (layer, start, end))
            assert layer == number and (start, end) == tuple(boundaries[number:number + 2])
            assert start < end <= count and end - start <= 64
            body, _ = function(source, name)
            tuples = re.findall(r'EvaluateOneStatic_0<(\d+)U, (\d+)U>\((\d+)U, (\d+)U, (\d+)U, (\d+)U,', body)
            assert len(tuples) == end - start
            for op, regs, destination, a, b, c in tuples:
                assert int(regs) == registers
                literal.extend(map(int, (op, destination, a, b, c)))
        assert literal == words[instruction_base:layer_base]
        prefix = expected[4]
        if prefix is not None:
            assert prefix in boundaries and 0 < prefix < count
        result.append(dict(offset=offset, count=count, registers=registers, prefix=prefix,
            roots=roots, boundaries=boundaries, layers=[row[3] for row in layers],
            program_end=layer_base + len(layers) + 1))
    for previous, following in zip(result, result[1:]):
        assert previous['program_end'] == following['offset']
    assert result[-1]['program_end'] == len(words)
    return result


def read_root(context, index, variable, indent):
    # This is exactly ReadValue_0's five ordered representation loads. The
    # pinned ScratchAddress_0 is index*5+component; rowbase/prefix are unused.
    out = [indent + f'RetainedTriple_0 {variable};']
    for component, field in enumerate(('hi_0', 'lo_0', 'tail_0', 'error_0', 'valid_0')):
        word = f'{context}->retainedScratch_0[{index * 5 + component}U]'
        value = word if component == 4 else f'RPFromBits_0({word})'
        out.append(indent + f'{variable}.{field} = {value};')
    return out


def specialized_body(stage, programs):
    p = PROFILES[stage]
    context, lane, rowbase = p['context'], p['lane'], p['rowbase']
    lines = [BEGIN.rstrip(), '    if (lanemask() != 0U)', '    {',
        f'        const uint32_t SiriusCanonicalBase = 1U + LOAD_WORD({context}->globalParams_0->inputs_0, 0U) * {p["stride"]}U;']
    counts = dict(branches=0, layer_calls=0, validation_roots=0, publication_roots=0)
    for program in programs:
        requirements = (20, 100) if stage == 'endpoint' else (p['roots'],)
        for required in requirements:
            condition = [f'{p["base"]} == SiriusCanonicalBase + {program["offset"]}U',
                f'{p["count"]} == {program["count"]}U', f'{p["registers"]} == {program["registers"]}U', f'{lane} == 0U']
            if p['layers']:
                condition.append(f'{p["layers"]} == {len(program["layers"])}U')
            if stage != 'transport':
                condition.append(f'required_0 == {required}U')
            if stage == 'endpoint':
                condition.append(f'prefix_4 == {program["prefix"]}U')
            # Uniform all-of-active gating preserves the original mixed-program
            # mask schedule by leaving every such invocation in the fallback.
            lines += ['        if (all(' + ' && '.join(condition) + '))', '        {',
                f'            if ({lane} == 0U) {context}->retainedStatus_0 = 0U;']
            first = program['prefix'] if stage == 'endpoint' and required == 100 else 0
            needed = program['prefix'] if stage == 'endpoint' and required == 20 else program['count']
            assert first in program['boundaries'] and needed in program['boundaries']
            for number, name in enumerate(program['layers']):
                start, end = program['boundaries'][number:number + 2]
                if not first <= start < needed:
                    continue
                assert end <= needed
                lines += [f'            {name}({p["row"]}, {rowbase}, {context});',
                    f'            if ({context}->retainedStatus_0 != 0U) return false;']
                counts['layer_calls'] += 1
            # Preserve complete ordered validation before the first root write.
            for index in program['roots'][:required]:
                lines += ['            {', f'                if ({index}U >= {p["registers"]})', '                {']
                if stage == 'transport':
                    lines += [f'                    RejectRetainedProgram_0({context});', '                }', '                else', '                {',
                        f'                    if ({context}->retainedScratch_0[{index * 5 + 4}U] == 0U) RejectRetainedProgram_0({context});']
                else:
                    lines += ['                    return false;', '                }']
                    lines += read_root(context, index, 'SiriusRoot', '                ')
                    lines += ['                if (SiriusRoot.valid_0 == 0U) return false;']
                if stage == 'transport':
                    lines += ['                }']
                lines += ['            }']
                counts['validation_roots'] += 1
            if stage == 'transport':
                lines += [f'            if ({context}->retainedStatus_0 != 0U) return false;']
            lines += [f'            if ({lane} == 0U)', '            {']
            for number, index in enumerate(program['roots'][:required]):
                lines += ['                {'] + read_root(context, index, 'SiriusRoot', '                    ')
                address = f'{rowbase} + 1004U + stage_2 * 200U + {number * 5}U' if stage == 'transport' else f'{rowbase} + {4 + number * 5}U'
                lines += [f'                    WriteStored_0({address}, &SiriusRoot, {context});', '                }']
                counts['publication_roots'] += 1
            lines += ['            }', '            return true;', '        }']
            counts['branches'] += 1
    lines += ['    }', END.rstrip()]
    return '\n'.join(lines) + '\n', counts


def derive(stage):
    source_path = R / f'bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/{stage}.ispc'
    source = source_path.read_text()
    scratch, _ = function(source, 'ScratchAddress_0')
    assert re.search(r'return index_\d+ \* 5U \+ component_\d+;', scratch)
    programs = canonical_programs(source, stage)
    insertion, counts = specialized_body(stage, programs)
    original_function, opening = function(source, 'EvaluateProgram_0')
    assert BEGIN not in source and END not in source
    candidate = source[:opening] + '\n' + insertion + source[opening:]
    before = f'#define sirius_simd_{stage}_packet sirius_simd_v4_{stage}_packet'
    after = f'#define sirius_simd_{stage}_packet sirius_simd_canonical_{stage}_packet'
    assert candidate.count(before) == 1 and after not in candidate
    candidate = candidate.replace(before, after, 1)
    inverse = candidate.replace(after, before, 1).replace('\n' + insertion, '', 1)
    assert inverse == source
    assert candidate.count(original_function) == 0  # Opening alone gained the block.
    assert len(re.findall(r'EvaluateOneStatic_0<\d+U, \d+U>\(', candidate)) == sum(row['count'] for row in programs)
    return candidate, dict(stage=stage, source_inverse_exact=True, original_fallback_body_unchanged=True,
        opcode_layer_definitions_and_tuple_order_unchanged=True, **counts,
        canonical_programs=[{key: row[key] for key in ('offset', 'count', 'registers', 'prefix')} for row in programs])


def main():
    assert sys.argv[1:] == ['derive'], 'Only the fixed source derivation is supported; no force or compiler mode.'
    contract = json.loads((W / 'canonical_metadata_contract.json').read_text())
    assert pin(Path(__file__)) == contract['generator']
    for name, expected in contract['source_inputs'].items():
        assert pin(R / name) == expected, name
    generated = {stage: derive(stage) for stage in STAGES}
    assert not B.exists(), 'Never overwrite an earlier candidate or output.'
    for stage, (source, record) in generated.items():
        assert pin_bytes(source.encode()) == contract['derived_sources'][stage]
        assert record == contract['source_findings'][stage]
    target = B / 'ispc-source'
    target.mkdir(parents=True)
    for stage, (source, _record) in generated.items():
        with (target / (stage + '.ispc')).open('x') as stream:
            stream.write(source)
    print(json.dumps(dict(prepared_only=True, compiler_run=False, numerical_qualification_transferred=False,
        generated_sources={stage: str((target / (stage + '.ispc')).relative_to(R)) for stage in STAGES},
        findings=contract['source_findings'], prospective_emit_commands=contract['prospective_emit_commands'])))


if __name__ == '__main__':
    main()
