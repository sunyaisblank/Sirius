"""Bind fresh Transport component evidence to the exact current R objects."""
from pathlib import Path
import hashlib
import json


def check(root):
    def doc(path):
        return json.loads(path.read_text())
    def pin(path):
        with path.open('rb') as stream:
            return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}
    def named(row):
        path = root / row['path']
        assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}
        return path
    work = root / 'out/native-cpu-accumulate-compatibility'
    for phase in ('link', 'run'):
        owner = doc(work / phase / 'owner.json')
        assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    build = doc(work / 'compile-completion.json')
    assert build['passed'] and build['diagnostic_CPP_compilations'] == 1 and build['diagnostic_links'] == 1
    assert build['additional_CPP_compilations_in_link_continuation'] == 0
    original = doc(work / 'compile/owner.json')
    assert original['status'] == 'terminal' and not original['passed'] and original['child_exit'] == 1
    assert original['source_unchanged'] and original['whole_inputs_unchanged']
    assert not original['cleanup_errors'] and not original['remaining_owned_processes']
    result_path = work / 'result.json'
    result = doc(result_path)
    assert result['passed'] and result['runtime_exit'] == 0
    comparison = result['comparison']
    for key, expected in (('case_count', 38), ('active_rows', 38), ('capacity_rows', 40),
                          ('packets_per_implementation', 5), ('raw_calls', 10),
                          ('compared_output_words', 187960)):
        assert comparison[key] == expected
    assert comparison['input_unchanged'] and comparison['inactive_sentinels_preserved']
    assert comparison['full_word_match'] and comparison['first_mismatch'] is None
    baseline, candidate = result['output_hashes']['baseline'], result['output_hashes']['candidate']
    assert baseline == candidate and baseline['bytes'] == 751840
    provider = root / 'out/native-cpu-accumulate-inline/provider'
    contract = doc(provider / 'source-contract.json')
    compiled_path = root / 'out/native-cpu-accumulate-inline/ispc-inline-compile-result.json'
    compiled = doc(compiled_path)
    assert compiled['all_objects_compiled'] and compiled['all_prepared_source_bytes_unchanged']
    objects = contract['native_objects']
    assert len(objects) == 8 and result['all8_native_objects'] == objects
    old = doc(root / 'out/native-cpu-transport-worker-occupancy/provider/source-contract.json')
    assert result['original_baseline_object'] == old['native_objects'][4]
    transport, = [row['object'] for row in compiled['objects'] if row['stage'] == 'transport']
    assert result['candidate_transport_object'] == transport == objects[4]
    link_path, build_path, executable = (named(result[key]) for key in ('link_proof', 'compile_result', 'executable'))
    link = doc(link_path)
    assert link['passed'] and link['executable'] == result['executable']
    for item in link['objects']:
        assert item['all_retained_function_opcodes_match']
    assert len(link['objects']) == 2
    assert {tuple(sorted(item['object'].items())) for item in link['objects']} == {
        tuple(sorted(row.items())) for row in (transport, result['renamed_baseline_object'])}
    current = doc(provider / 'linked-native-proof.json')
    assert current['passed'] and [item['object'] for item in current['objects']] == objects
    assert all(item['all_retained_function_opcodes_match'] for item in current['objects'])
    for row in objects + [result['original_baseline_object'], result['renamed_baseline_object']]:
        named(row)
    inputs = {str(path.relative_to(root)): pin(path) for path in work.rglob('*') if path.is_file()}
    for path in (provider / 'linked-native-proof.json', provider / 'source-contract.json', compiled_path, executable):
        inputs[str(path.relative_to(root))] = pin(path)
    return {'inputs': inputs, 'scope': {'component_result': pin(result_path),
        'component_link_proof': pin(link_path), 'component_compile_result': pin(build_path),
        'current_R_link_proof': pin(provider / 'linked-native-proof.json'),
        'candidate_transport_object': transport, 'compared_output_words': 187960,
        'scope': 'Fresh old-SIMD/new-R direct Transport compatibility on all38 frozen inputs/full capacity40, tied to actual standalone and provider R object code. Finite implementation compatibility only; no independent-math, all-domain, benefit, frame or scientific/release claim.'}}
