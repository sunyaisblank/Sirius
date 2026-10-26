"""Require the source-bound, completed Endpoint component before composition."""
from pathlib import Path
import hashlib
import json

R = Path.cwd()
D = R / 'out/native-cpu-endpoint-header-domain'
OLD = R / 'out/native-cpu-v4-provider/provider'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def checked(path, expected):
    assert pin(path) == {key: expected[key] for key in ('bytes', 'sha256')}, path


def terminal(path):
    owner = doc(path)
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0, path
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'], path
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors'], path
    return owner


def component():
    derivation = doc(D / 'derivation.json')
    for name, expected in derivation['basis'].items():
        checked(R / name, expected)
    for name, expected in derivation['derived'].items():
        checked(D / name, expected)
    assert set(derivation['derived']) == {'native_packet.hpp', 'admission.cpp', 'endpoint.cpp'}
    assert derivation['admitted_program_difference'] == dict(stage='Endpoint', relative_word=11241,
        bytes=4, canonical=612, allowed='Any uint32; the unchanged shader owns acceptance/refusal')
    terminal(D / 'build/owner.json')
    build = doc(D / 'build-result.json')
    assert build['new_compilations'] == 2 and len(build['objects']) == 2
    objects = [row['object'] for row in build['objects']]
    assert {row['path'] for row in objects} == {
        'bin/linux-gcc/native-cpu-endpoint-header-domain/endpoint.o',
        'bin/linux-gcc/native-cpu-endpoint-header-domain/admission.o'}
    for record in build['objects']:
        assert record['exit'] == 0
        checked(R / record['object']['path'], record['object'])
    proof = doc(D / 'code-proof.json')
    assert proof['passed'] and len(proof['objects']) == 2
    assert sum(row['unchanged_function_count'] for row in proof['objects']) == 168
    old = doc(OLD / 'source-contract.json')
    candidates = {Path(row['path']).name: row for row in objects}
    for row, index, name, symbol, count in zip(proof['objects'], (1, 3),
            ('endpoint.o', 'admission.o'),
            ('sirius_native_endpoint_try_packet', 'sirius_native_endpoint_admit_packet'), (166, 2)):
        assert row['original'] == old['native_objects'][index]
        assert row['candidate'] == candidates[name]
        assert row['changed_functions'] == [symbol] and row['unchanged_function_count'] == count
        assert row['function_bytes_and_relocation_targets_identical']
        assert row['optional_ISA_instructions'] == 0
        checked(R / row['original']['path'], row['original'])
    terminal(D / 'compatibility/owner.json')
    compatibility = doc(D / 'compatibility.json')
    assert compatibility['passed'] and compatibility['inputs_unchanged'] and compatibility['source_unchanged']
    assert compatibility['objects'] == objects
    comparison = compatibility['comparison']
    assert comparison['passed']
    for key, value in dict(selection_rows=23, capacity=32, header_values=6, layouts=2).items():
        assert comparison[key] == value, key
    for key in ('whole3829_words_and_inactive_capacity_and_canaries', 'input_immutable',
                'MXCSR_controls_unchanged', 'caller_MXCSR_state_restored'):
        assert comparison[key] is True, key
    for key in ('reference_rows', 'scalar_packets', 'raw_packets', 'pure_calls',
                'other_program_word_refusals', 'boundary_refusals',
                'compared_output_words', 'compared_input_words'):
        assert type(comparison[key]) is int and comparison[key] > 0, key
    return objects


def contract_inputs(contract):
    for name, expected in contract['source_inputs'].items():
        checked(R / name, expected)
