"""Require this actual field Transport's complete fixture and exact code joins."""
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
        assert pin(path) == {k: row[k] for k in ('bytes', 'sha256')}
        return path
    work = root / 'out/native-cpu-field-compatibility'
    for phase in ('proof', 'run'):
        owner = doc(work / phase / 'owner.json')
        assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
        assert owner['whole_inputs_unchanged'] and owner['source_unchanged']
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    result = doc(work / 'result.json'); comparison = result['comparison']
    assert result['passed'] and result['runtime_exit'] == 0 and result['all_prepared_inputs_unchanged']
    for key, value in (('case_count',38), ('active_rows',38), ('capacity_rows',40),
                       ('packets_per_implementation',5), ('raw_calls',10), ('compared_output_words',187960)):
        assert comparison[key] == value
    assert comparison['full_word_match'] and comparison['inactive_sentinels_preserved'] and comparison['input_unchanged']
    assert comparison['mismatch_count'] == 0 and comparison['first_mismatch'] is None
    assert result['output_hashes']['baseline'] == result['output_hashes']['candidate']
    provider = root / 'out/native-cpu-field-simd/provider'
    contract = doc(provider / 'source-contract.json'); current = doc(provider / 'linked-native-proof.json')
    objects = contract['native_objects']; assert len(objects) == 11
    assert current['passed'] and [row['object'] for row in current['objects']] == objects
    assert all(row['all_retained_function_opcodes_match'] for row in current['objects'])
    assert result['transport_component_objects'] == [objects[4],objects[7],objects[3]]
    assert result['candidate_transport_object'] == objects[7]
    build = doc(named(result['compile_result'])); link = doc(named(result['link_proof']))
    assert build['passed'] and build['diagnostic_shim_compilations'] == 1 and build['new_links'] == 1
    assert build['additional_compilations_in_proof_completion'] == build['additional_links_in_proof_completion'] == 0
    assert link['passed'] and len(link['objects']) == 4
    assert [row['object'] for row in link['objects']] == [result['renamed_baseline_object'],objects[4],objects[7],build['shim_object']]
    assert all(row['all_retained_function_opcodes_match'] for row in link['objects'])
    for row in objects + [result['renamed_baseline_object'], build['driver_object'], build['shim_object']]:
        named(row)
    executable = named(result['executable'])
    inputs = {str(path.relative_to(root)): pin(path) for path in work.rglob('*') if path.is_file()}
    for path in (provider / 'source-contract.json', provider / 'linked-native-proof.json', executable):
        inputs[str(path.relative_to(root))] = pin(path)
    return {'inputs': inputs, 'scope': {'component_result': pin(work / 'result.json'),
        'component_code_proof': result['link_proof'], 'current_provider_code_proof': pin(provider / 'linked-native-proof.json'),
        'candidate_transport_objects': [objects[4],objects[7]], 'compared_output_words':187960,
        'scope':'Actual full38 field Transport, unchanged fixture driver/baseline, all187960words. Exact standalone and current provider code joins. Finite compiler/executor compatibility; no independent math or full scientific/release claim.'}}
