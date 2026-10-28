"""Closed private Transport/Endpoint metadata with original mixed-lane order."""
from pathlib import Path
import hashlib
import json
import re
import sys
import types

R = Path(__file__).resolve().parents[2]
V = Path(__file__).resolve().parent
STAGES = ('transport', 'endpoint')
METADATA = V / 'canonical_metadata.py'
METADATA_SHA = '4136ab9b4a666c86d1a7ff7837c7f6d89c7036ba185163dfeb4f262f9410ef23'
PINS = {
    'transport': {
        'AS': dict(bytes=626976, sha256='0fc91ff0b756fddd055d5c0322f5eece0557a64799ef31893933be33cdba10f6'),
        'publication': dict(bytes=771952, sha256='761703bd80a16343ac894cec29a11a0ff2211d550f444e48b872846330011a11')},
    'endpoint': {
        'AS': dict(bytes=795302, sha256='d8ae7cdde50211537678ec197f18c374edd2cec653974d9d0b7c33c50db8c1e1'),
        'publication': dict(bytes=1242957, sha256='5728c1ecff0723ac68e97c13c82269c024a5b15cbb3870c29c92fae7c770c4c8')},
}


def pin_bytes(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def authority():
    assert pin(METADATA)['sha256'] == METADATA_SHA
    module = types.ModuleType('pinned_canonical_metadata')
    module.__file__ = str(METADATA)
    exec(compile(METADATA.read_text(), str(METADATA), 'exec'), module.__dict__)
    return module


def target(stage):
    return R / f'bin/linux-gcc/native-cpu-canonical-closed-{stage}/ispc-source/{stage}.ispc'


def eligibility(stage, program, layer):
    if layer >= len(program['layers']):
        return 'false'
    family = 'SiriusFamily' + str(program['family'])
    if stage == 'transport':
        return family
    start, end = program['boundaries'][layer:layer + 2]
    prefix = program['prefix']
    assert end <= prefix or start >= prefix
    return family + (' && required_0 == 20U' if end <= prefix else ' && required_0 == 100U')


def validations(stage, programs, p):
    ctx, regs = p['context'], p['registers']
    count = p['roots']
    out = []
    for number in range(count):
        outer = stage == 'endpoint' and number >= 20
        out += ['    if (required_0 == 100U)' if outer else '    /* Both Endpoint phases share this root position. */' if stage == 'endpoint' else '    /* Original Transport root position. */', '    {',
            '        bool SiriusRootInvalid = false;']
        for family, program in enumerate(programs):
            root = program['roots'][number]
            out += [f'        if (SiriusFamily{family})', '        {',
                f'            if ({root}U >= {regs})', '            {',
                '                SiriusRootInvalid = true;', '            }', '            else', '            {']
            if stage == 'transport':
                out += [f'                SiriusRootInvalid = {ctx}->retainedScratch_0[{root * 5 + 4}U] == 0U;']
            else:
                # Exact original ReadValue representation loads, no address escape.
                out += ['                RetainedTriple_0 SiriusRoot;']
                for component, field in enumerate(('hi_0', 'lo_0', 'tail_0', 'error_0', 'valid_0')):
                    word = f'{ctx}->retainedScratch_0[{root * 5 + component}U]'
                    value = word if component == 4 else f'RPFromBits_0({word})'
                    out += [f'                SiriusRoot.{field} = {value};']
                out += ['                SiriusRootInvalid = SiriusRoot.valid_0 == 0U;']
            out += ['            }', '        }']
        if stage == 'transport':
            # One original union-of-invalid-lanes status update per root position.
            out += [f'        if (SiriusRootInvalid) RejectRetainedProgram_0({ctx});']
        else:
            out += ['        if (SiriusRootInvalid) return false;']
        out += ['    }']
    if stage == 'transport':
        out += [f'    if ({ctx}->retainedStatus_0 != 0U) return false;']
    return out


def publications(stage, programs, p):
    ctx, rowbase = p['context'], p['rowbase']
    out = ['    /* All required roots were validated before this first store. */',
        f'    if ({p["lane"]} == 0U)', '    {']
    for number in range(p['roots']):
        outer = stage == 'endpoint' and number >= 20
        if outer:
            out += ['        if (required_0 == 100U)']
        out += ['        {']
        for component in range(5):
            out += [f'            uint32_t SiriusPublishedWord{component} = 0U;']
        for family, program in enumerate(programs):
            root = program['roots'][number]
            out += [f'            if (SiriusFamily{family})', '            {']
            for component in range(5):
                out += [f'                SiriusPublishedWord{component} = {ctx}->retainedScratch_0[{root * 5 + component}U];']
            out += ['            }']
        address = f'{rowbase} + 1004U + stage_2 * 200U + {number * 5}U' if stage == 'transport' else f'{rowbase} + {4 + number * 5}U'
        # Preserve original word-major combined-lane publication: five reads
        # from each active family first, then five stores under the union mask.
        for component in range(5):
            field = address + (f' + {component}U' if component else '')
            out += [f'            *(&OUTPUT_WORD({ctx}, {field})) = SiriusPublishedWord{component};']
        out += ['        }']
    out += ['    }']
    return out


def closed_function(stage, original, programs, p):
    signature = original.split('\n{', 1)[0]
    ctx, lane = p['context'], p['lane']
    out = [signature, '{', '    /* Closed positive-job private path; mixed families/phases remain varying. */',
        '    if (lanemask() == 0U) return false;',
        f'    const uint32_t SiriusCanonicalBase = 1U + LOAD_WORD({ctx}->globalParams_0->inputs_0, 0U) * {p["stride"]}U;']
    for family, program in enumerate(programs):
        condition = [f'{p["base"]} == SiriusCanonicalBase + {program["offset"]}U',
            f'{p["count"]} == {program["count"]}U', f'{p["registers"]} == {program["registers"]}U',
            f'{p["layers"]} == {len(program["layers"])}U']
        if stage == 'endpoint':
            condition.append(f'prefix_4 == {program["prefix"]}U')
        out += [f'    bool SiriusFamily{family} = ' + ' && '.join(condition) + ';']
    bad = f'!(SiriusFamily0 || SiriusFamily1) || {lane} != 0U'
    if stage == 'endpoint':
        bad += ' || (required_0 != 20U && required_0 != 100U)'
    out += [f'    if ({bad}) return false;', f'    if ({lane} == 0U) {ctx}->retainedStatus_0 = 0U;']
    order = []
    for layer in range(max(len(program['layers']) for program in programs)):
        a, b = [eligibility(stage, program, layer) for program in programs]
        out += [f'    /* Original lockstep layer {layer}; family0 then family1, then one status check. */', '    {',
            f'        bool SiriusExecute0 = {a};', f'        bool SiriusExecute1 = {b};']
        for family, program in enumerate(programs):
            if layer >= len(program['layers']):
                continue
            name = program['layers'][layer]
            out += [f'        if (SiriusExecute{family})', '        {',
                f'            {name}({p["row"]}, {p["rowbase"]}, {ctx});', '        }']
            order.append([layer, family, name])
        out += ['        if (SiriusExecute0 || SiriusExecute1)', '        {',
            f'            if ({ctx}->retainedStatus_0 != 0U) return false;', '        }', '    }']
    out += validations(stage, programs, p) + publications(stage, programs, p)
    out += ['    return true;', '}']
    text = '\n'.join(out)
    assert 'all(' not in text and 'for (' not in text and 'for(' not in text
    assert text.count('LOAD_WORD(') == 1
    assert not any(name in text for name in ('EvaluateKnownLayer_0', 'EvaluateOne_0(', 'ReadValue_0(', 'WriteStored_0('))
    assert text.rindex('if (SiriusRootInvalid)') < text.index('*(&OUTPUT_WORD(')
    assert text.count('*(&OUTPUT_WORD(') == p['roots'] * 5
    return text, order


def caller_proof(stage, source, metadata):
    assert source.count('EvaluateProgram_0(') == source.count('ComputeOwnedRow_0(') == 2
    compute, _ = metadata.function(source, 'ComputeOwnedRow_0')
    export = source[source.index(f'\nexport void sirius_simd_{stage}_packet('):]
    assert re.search(r'for \(uniform uint32_t invocation = 0; invocation < 64U; \+\+invocation\) \{\s*'
        r'if \(invocation == 0U\) \{\s*varying KernelContext_0 context;\s*'
        r'memset\(\(void \* uniform\)&context, 0, sizeof\(context\)\);\s*'
        r'context.globalParams_0 = &params;\s*'
        r'RowID id; id.x = first_x \+ programIndex; id.y = 0U; id.z = 0U;\s*'
        r'if \(programIndex < groups\) ComputeOwnedRow_0\(id, invocation, &context\);\s*\}\s*\}', export)
    if stage == 'transport':
        expected = 'EvaluateProgram_0(row_4, rowBase_8, programBase_3, count_2, registers_2, layers_1, stage_3, lane_2, kernelContext_13)'
        assert compute.count(expected) == 1
        for line in ('layers_1 = 42U;', 'layers_1 = 117U;', 'programBase_3 = programBase_4;',
            'programBase_3 = programBase_2;', 'count_2 = count_3;', 'count_2 = count_1;'):
            assert compute.count(line) == 1
    else:
        expected = 'EvaluateProgram_0(row_2, rowBase_11, _S402, _S403, registers_3, (uint32_t)(_S404), lane_1, layers_2, prefix_5, kernelContext_13)'
        assert compute.count(expected) == 1
        cleaned = '\n'.join(line for line in compute.splitlines() if not line.strip().startswith('#line'))
        assert re.search(r'if\(phase_index_0 == 0U\)\s*\{\s*_S404 = \(int\)\(20\);\s*\}\s*else\s*\{\s*_S404 = \(int\)\(100\);\s*\}', cleaned)
        select, _ = metadata.function(source, 'SelectEndpointProgram_0')
        for line in ('*layers_0 = 84U;', '*prefix_0 = 202U;',
            'if(LOAD_WORD(kernelContext_2->globalParams_0->inputs_0, *programBase_0 + 1U) == registers_0)'):
            assert select.count(line) == 1
        assert compute.count('uint32_t layers_2 = 111U;') == 1 and compute.count('uint32_t prefix_5 = 250U;') == 1


def derive(stage):
    assert stage in STAGES
    metadata = authority()
    original_path = R / f'bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/{stage}.ispc'
    publication_path = R / f'bin/linux-gcc/native-cpu-canonical-metadata/ispc-source/{stage}.publication.ispc'
    assert pin(original_path) == PINS[stage]['AS'] and pin(publication_path) == PINS[stage]['publication']
    original, publication = original_path.read_text(), publication_path.read_text()
    old_body, _ = metadata.function(original, 'EvaluateProgram_0')
    published_body, _ = metadata.function(publication, 'EvaluateProgram_0')
    as_binding = f'#define sirius_simd_{stage}_packet sirius_simd_v4_{stage}_packet'
    pub_binding = f'#define sirius_simd_{stage}_packet sirius_simd_publication_{stage}_packet'
    new_binding = f'#define sirius_simd_{stage}_packet sirius_simd_closed_{stage}_packet'
    assert publication.count(published_body) == publication.count(pub_binding) == original.count(old_body) == 1
    assert publication.replace(published_body, old_body, 1).replace(pub_binding, as_binding, 1) == original
    programs = metadata.canonical_programs(original, stage)
    assert len(programs) == 2
    for family, program in enumerate(programs):
        program['family'] = family
    p = metadata.PROFILES[stage]
    new_body, order = closed_function(stage, old_body, programs, p)
    candidate = publication.replace(published_body, new_body, 1).replace(pub_binding, new_binding, 1)
    assert candidate.replace(new_body, published_body, 1).replace(new_binding, pub_binding, 1) == publication
    assert candidate.replace(new_body, old_body, 1).replace(new_binding, as_binding, 1) == original
    assert new_body.count('EvaluateStaticLayer_') == sum(len(program['layers']) for program in programs)
    assert re.findall(r'\b(EvaluateStaticLayer_\d+_\d+_0)\(', new_body) == [row[2] for row in order]
    assert len(re.findall(r'EvaluateOneStatic_0<\d+U, \d+U>\(', candidate)) == sum(program['count'] for program in programs)
    caller_proof(stage, original, metadata)
    proof = dict(source_inverse_to_publication_exact=True, source_inverse_to_AS_exact=True,
        unchanged_outside_export_binding_and_EvaluateProgram=True, generic_fallback_removed_only_from_new_private_body=True,
        canonical_tuple_records_match_exact_tables=True, original_tuple_records=sum(program['count'] for program in programs),
        complete_context_memset_and_full8_physical_packing_unchanged=True, supported_callers_lane0=True,
        mixed_families_supported=True, homogeneous_all_active_guard_absent=True,
        per_layer_family_order=order, original_layer_call_sites=len(order),
        required_roots=40 if stage == 'transport' else [20, 100], root_positions=p['roots'],
        all_root_validation_before_any_publication=True, original_combined_mask_five_word_store_order=True,
        publication_word_store_sites=p['roots'] * 5,
        programs=[{key: program[key] for key in ('offset', 'count', 'registers', 'prefix', 'roots', 'boundaries')} for program in programs])
    if stage == 'endpoint':
        proof.update(live_Endpoint11241_guard_unchanged=True,
            Endpoint_prefix_state_and_whole_state_rollback_unchanged=True,
            required20_and100_supported_per_lane=True)
    return candidate, proof


def main():
    assert sys.argv[1:] == ['derive'], 'Only fixed source preparation; no force/compile/runtime mode.'
    contract = json.loads((V / 'canonical_private_remaining_contract.json').read_text())
    assert pin(Path(__file__)) == contract['generator']
    for name, expected in contract['source_inputs'].items():
        assert pin(R / name) == expected, name
    generated = {stage: derive(stage) for stage in STAGES}
    for stage, (candidate, proof) in generated.items():
        assert pin_bytes(candidate.encode()) == contract['derived_sources'][stage]
        assert proof == contract['source_proofs'][stage]
        assert not target(stage).parents[1].exists(), 'Do not overwrite an earlier closed candidate.'
    for stage, (candidate, _proof) in generated.items():
        path = target(stage)
        path.parent.mkdir(parents=True)
        with path.open('x') as stream:
            stream.write(candidate)
    print(json.dumps(dict(prepared_only=True, compiler_run=False, numerical_qualification_transferred=False,
        outputs={stage: dict(path=str(target(stage).relative_to(R)), **pin(target(stage))) for stage in STAGES},
        prospective_emit_commands=contract['prospective_emit_commands'])))


if __name__ == '__main__':
    main()
