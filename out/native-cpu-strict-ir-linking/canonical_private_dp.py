"""Prepare a closed private DP canonical source; never compile or execute it."""
from pathlib import Path
import hashlib
import json
import re
import sys

R = Path(__file__).resolve().parents[2]
W = Path(__file__).resolve().parent
INPUT = R / 'bin/linux-gcc/native-cpu-canonical-metadata/ispc-source/dopri_phase.publication.ispc'
INPUT_PIN = dict(bytes=756241, sha256='c4c1936fa1bb19e82a0f4ff388335ffb985b7ddfa8e576754bc5ad2574c68b75')
B = R / 'bin/linux-gcc/native-cpu-canonical-closed-dp'
OUTPUT = B / 'ispc-source/dopri_phase.ispc'
BEGIN = '    /* BEGIN pinned canonical metadata specialization. */\n'
END = '    /* END pinned canonical metadata specialization; exact fallback follows. */\n'
SIGNATURE = ('static bool EvaluateProgram_0(uint32_t row_2, uniform uint32_t rowBase_6, '
    'uint32_t programBase_1, uint32_t count_0, uint32_t registers_1, uint32_t required_0, '
    'uniform uint32_t lane_0, varying KernelContext_0 * uniform kernelContext_8)\n')
BASE = 'const uint32_t SiriusCanonicalBase = 1U + LOAD_WORD(kernelContext_8->globalParams_0->inputs_0, 0U) * 1810U;'
PREDICATE = ('programBase_1 == SiriusCanonicalBase + 0U && count_0 == 2088U && '
    'registers_1 == 629U && lane_0 == 0U && required_0 == 200U')
BEFORE = '#define sirius_simd_dopri_phase_packet sirius_simd_publication_dopri_phase_packet'
AFTER = '#define sirius_simd_dopri_phase_packet sirius_simd_closed_dp_packet'


def pin_bytes(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def derive():
    assert pin(INPUT) == INPUT_PIN
    source = INPUT.read_text()
    assert source.count(SIGNATURE) == source.count(BEGIN) == source.count(END) == 1
    start = source.index(SIGNATURE)
    finish = source.index('\n}\n\nstatic void ComputeOwnedRow_0(', start) + 2
    original = source[start:finish]
    prefix = (SIGNATURE + '{\n' + BEGIN + '    if (lanemask() != 0U)\n    {\n'
        + '        ' + BASE + '\n        if (all(' + PREDICATE + '))\n        {\n')
    suffix = '        }\n    }\n' + END
    assert original.startswith(prefix)
    fast_end = original.index(END) + len(END)
    fast = original[:fast_end]
    assert fast.endswith(suffix)
    payload = fast[len(prefix):-len(suffix)]
    assert payload.endswith('            return true;\n')
    assert all(not line or line.startswith('        ') for line in payload.splitlines())
    payload = '\n'.join(line[8:] for line in payload.splitlines()) + '\n'
    # Only unwrap the already checked branch. No payload/tuple/FP rewrite.
    new_function = (SIGNATURE + '{\n'
        + '    /* Private closed DP: same positive immutable-program job and original FP scope required. */\n'
        + '    if (lanemask() == 0U) return false;\n    ' + BASE + '\n'
        + '    if (!all(' + PREDICATE + ')) return false;\n' + payload + '}')
    assert new_function.count('EvaluateStaticLayer_0_') == 45
    calls = re.findall(r'\bEvaluateStaticLayer_0_(\d+)_0\(row_2, rowBase_6, kernelContext_8\);', payload)
    assert list(map(int, calls)) == list(range(45))
    assert new_function.count('if (kernelContext_8->retainedStatus_0 != 0U) return false;') == 45
    assert new_function.count('if (SiriusRoot.valid_0 == 0U) return false;') == 200
    assert new_function.count('const uint32_t SiriusPublishedWord0 = ') == 200
    assert new_function.count('*(&OUTPUT_WORD(') == 1000
    assert not any(word in new_function for word in (
        'EvaluateKnownLayer_0', 'EvaluateOne_0(', 'ReadValue_0(', 'WriteStored_0(', 'for(', 'for (', 'while ('))
    # No publication can precede the last root validation.
    assert new_function.rindex('if (SiriusRoot.valid_0 == 0U) return false;') < new_function.index('*(&OUTPUT_WORD(')
    assert new_function.count('LOAD_WORD(') == 1  # Only the immutable capacity used for the guard base.
    assert payload.count('retainedScratch_0[') == 2000
    # Verify the root list/addresses independently against the exact admitted table.
    table = R / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch/dopri_phase.program.hpp'
    match = re.search(r'dopri_phase_program\[\] = \{(.*?)\};', table.read_text(), re.S)
    assert match and not re.sub(r'\d+u|[\s,]', '', match[1])
    words = [int(value) for value in re.findall(r'(\d+)u', match[1])]
    assert len(words) == 10689 and words[:3] == [2088, 629, 200]
    tuples = re.findall(r'EvaluateOneStatic_0<(\d+)U, (\d+)U>\((\d+)U, (\d+)U, (\d+)U, (\d+)U,', source)
    assert len(tuples) == 2088
    literal = []
    for op, registers, destination, a, b, c in tuples:
        assert int(registers) == 629
        literal.extend(map(int, (op, destination, a, b, c)))
    assert literal == words[203:10643]
    roots = words[3:203]
    for number, root in enumerate(roots):
        assert root < 629
        validation = (f'SiriusRoot.valid_0 = kernelContext_8->retainedScratch_0[{root * 5 + 4}U];\n'
            '        if (SiriusRoot.valid_0 == 0U) return false;')
        assert validation in payload
        publication = [f'const uint32_t SiriusPublishedWord{component} = kernelContext_8->retainedScratch_0[{root * 5 + component}U];'
            for component in range(5)]
        cursor = 0
        for line in publication:
            cursor = payload.index(line, cursor) + len(line)
        for component in range(5):
            address = f'rowBase_6 + {4 + number * 5}U' + (f' + {component}U' if component else '')
            line = f'*(&OUTPUT_WORD(kernelContext_8, {address})) = SiriusPublishedWord{component};'
            cursor = payload.index(line, cursor) + len(line)
    # Supported callers are the original owned-row export, invocation0 only.
    call = 'EvaluateProgram_0(row_3, rowBase_7, programBase_2, count_1, registers_2, 200U, lane_1, kernelContext_9)'
    assert source.count('EvaluateProgram_0(') == 2 and source.count(call) == 1
    assert source.count('ComputeOwnedRow_0(') == 2
    compute = source[finish:source.index('\nexport void sirius_simd_dopri_phase_packet(', finish)]
    for line in (
        'uint32_t programBase_2 = 1U + 1810U * rows_0;',
        'uint32_t count_1 = LOAD_WORD(kernelContext_9->globalParams_0->inputs_0, programBase_2);',
        'uint32_t registers_2 = LOAD_WORD(kernelContext_9->globalParams_0->inputs_0, programBase_2 + 1U);'):
        assert compute.count(line) == 1
    export = source[source.index('\nexport void sirius_simd_dopri_phase_packet(', finish):]
    assert re.search(r'for \(uniform uint32_t invocation = 0; invocation < 64U; \+\+invocation\) \{\s*'
        r'if \(invocation == 0U\) \{\s*varying KernelContext_0 context;\s*'
        r'memset\(\(void \* uniform\)&context, 0, sizeof\(context\)\);\s*'
        r'context.globalParams_0 = &params;\s*'
        r'RowID id; id.x = first_x \+ programIndex; id.y = 0U; id.z = 0U;\s*'
        r'if \(programIndex < groups\) ComputeOwnedRow_0\(id, invocation, &context\);\s*\}\s*\}', export)
    candidate = source[:start] + new_function + source[finish:]
    assert candidate.count(BEFORE) == 1 and AFTER not in candidate
    candidate = candidate.replace(BEFORE, AFTER, 1)
    restored = candidate.replace(AFTER, BEFORE, 1)
    assert restored.count(new_function) == 1
    restored = restored.replace(new_function, original, 1)
    assert restored == source
    assert candidate.count('EvaluateOneStatic_0<') == source.count('EvaluateOneStatic_0<')
    assert len(re.findall(r'EvaluateOneStatic_0<\d+U, \d+U>\(', candidate)) == 2088
    return candidate, dict(source_inverse_exact=True, unchanged_outside_export_binding_and_EvaluateProgram=True,
        old_generic_body=pin_bytes(original.encode()), new_closed_body=pin_bytes(new_function.encode()),
        canonical_payload_unwrapped_only=True, layer_calls=45, ordered_layer_numbers=list(range(45)),
        layer_status_checks=45, validated_roots=200, published_roots=200, publication_words=1000,
        all_validation_before_first_publication=True, complete_private_context_memset_preserved=True,
        original_owned_row_shape_IO_guards_and_callers_unchanged=True, original_tuple_records=2088,
        canonical_tuple_records_match_exact_admitted_table=True, canonical_root_validation_and_publication_match_table=True,
        sole_supported_caller_required200_lane0=True, full8_physical_payloads_and_tags_unchanged=True,
        generic_program_fallback_removed_from_new_private_body_only=True)


def main():
    assert sys.argv[1:] == ['derive'], 'Only fixed source derivation; no force/compiler/runtime mode.'
    contract = json.loads((W / 'canonical_private_dp_contract.json').read_text())
    assert pin(Path(__file__)) == contract['generator']
    for name, expected in contract['source_inputs'].items():
        assert pin(R / name) == expected, name
    candidate, proof = derive()
    assert pin_bytes(candidate.encode()) == contract['derived_source'] and proof == contract['source_proof']
    assert not B.exists(), 'Do not overwrite an earlier candidate.'
    OUTPUT.parent.mkdir(parents=True)
    with OUTPUT.open('x') as stream:
        stream.write(candidate)
    print(json.dumps(dict(prepared_only=True, compiler_run=False, numerical_qualification_transferred=False,
        output=dict(path=str(OUTPUT.relative_to(R)), **pin(OUTPUT)), source_proof=proof,
        prospective_emit_command=contract['prospective_emit_command'])))


if __name__ == '__main__':
    main()
