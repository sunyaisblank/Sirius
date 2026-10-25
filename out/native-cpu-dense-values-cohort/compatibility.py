"""Require fresh full38 W execution and exact standalone/provider code associations."""
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
    work = root / 'out/native-cpu-dense-values-compatibility'
    for phase in ('link', 'run'):
        owner = doc(work / phase / 'owner.json')
        assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    result_path = work / 'result.json'
    result = doc(result_path)
    assert result['passed'] and result['runtime_exit'] == 0 and result['all_prepared_inputs_unchanged']
    comparison = result['comparison']
    for key, expected in (('case_count', 38), ('active_rows', 38), ('capacity_rows', 40),
                          ('packets_per_implementation', 5), ('raw_calls', 10), ('compared_output_words', 187960)):
        assert comparison[key] == expected
    assert comparison['input_unchanged'] and comparison['inactive_sentinels_preserved']
    assert comparison['full_word_match'] and comparison['first_mismatch'] is None
    assert result['output_hashes']['baseline'] == result['output_hashes']['candidate']
    assert result['output_hashes']['candidate']['bytes'] == 751840
    provider = root / 'out/native-cpu-dense-values/provider'
    contract = doc(provider / 'source-contract.json')
    compiled_path = root / 'out/native-cpu-dense-values/ispc-values-compile-result.json'
    compiled = doc(compiled_path)
    assert compiled['all_objects_compiled'] and compiled['all_prepared_source_bytes_unchanged']
    objects = contract['native_objects']
    assert len(objects) == 8 and result['all8_native_objects'] == objects
    baseline = doc(root / 'out/native-cpu-transport-worker-occupancy/provider/source-contract.json')['native_objects'][4]
    assert result['original_baseline_object'] == baseline
    transport, = [row['object'] for row in compiled['objects'] if row['stage'] == 'transport']
    assert result['candidate_transport_object'] == transport == objects[4]
    link_path, build_path, executable = (named(result[key]) for key in ('link_proof', 'compile_result', 'executable'))
    build = doc(build_path)
    assert build['passed'] and build['new_CPP_compilations'] == 0 and build['new_links'] == 1
    reused = doc(named(build['reused_R_build']))
    assert reused['passed'] and reused['driver_object'] == build['driver_object']
    assert reused['renamed_baseline_object'] == build['renamed_baseline_object']
    link = doc(link_path)
    assert link['passed'] and link['executable'] == result['executable'] and len(link['objects']) == 2
    assert all(item['all_retained_function_opcodes_match'] for item in link['objects'])
    assert {tuple(sorted(item['object'].items())) for item in link['objects']} == {
        tuple(sorted(row.items())) for row in (transport, result['renamed_baseline_object'])}
    current = doc(provider / 'linked-native-proof.json')
    assert current['passed'] and [item['object'] for item in current['objects']] == objects
    assert all(item['all_retained_function_opcodes_match'] for item in current['objects'])
    for row in objects + [baseline, result['renamed_baseline_object'], build['driver_object']]:
        named(row)
    inputs = {str(path.relative_to(root)): pin(path) for path in work.rglob('*') if path.is_file()}
    for path in (provider / 'linked-native-proof.json', provider / 'source-contract.json', compiled_path, executable,
                 named(build['driver_object']), named(build['reused_R_build'])):
        inputs[str(path.relative_to(root))] = pin(path)
    return {'inputs': inputs, 'scope': {'component_result': pin(result_path),
        'component_link_proof': pin(link_path), 'component_compile_result': pin(build_path),
        'current_AD_link_proof': pin(provider / 'linked-native-proof.json'),
        'candidate_transport_object': transport, 'compared_output_words': 187960,
        'scope': 'Fresh old/new AA raw Transport compatibility, all38 unchanged inputs/all187960words. Reused exact already-compiled R driver and renamed Q baseline; one new link/one workload, no CPP or old-math recompilation. Exact standalone/provider AA object-code association; finite compatibility only, no independent-math, benefit, frame or science/release completion claim.'}}
