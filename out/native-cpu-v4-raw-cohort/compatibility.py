"""Reuse the sealed genuine R38/Endpoint checks and join their exact AQ parts.

No unchanged fixture is relinked or rerun. Current provider code composition and
the fresh C8 selection/actual DopriPhase execution remain separate required gates.
"""
from pathlib import Path
import hashlib
import json
import zipfile


def check(root):
    def doc(path):
        return json.loads(path.read_text())
    def pin(path):
        with path.open('rb') as stream:
            return dict(bytes=path.stat().st_size,
                        sha256=hashlib.file_digest(stream, 'sha256').hexdigest())
    manifest_path = root / 'attestations/qualification/984ba457/native-cpu-v4-raw-packets-finite-compatibility/MANIFEST.json'
    assert pin(manifest_path)['sha256'] == 'b41bc17d44c4757c51415a1282da22692c85cbd589e92c11724e47da01cdc19b'
    manifest = doc(manifest_path)
    assert manifest['finite_component_passed'] and not manifest['provider_adopted']
    archive = root / manifest['archive']['path']
    assert pin(archive) == {key: manifest['archive'][key] for key in ('bytes', 'sha256')}
    with zipfile.ZipFile(archive) as zipped:
        def member(name):
            key = 'out/native-cpu-v4-raw-packets/' + name
            data = zipped.read(key)
            assert dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest()) == manifest['members'][key]
            return json.loads(data)
        finite = member('finite-result.json')
        code = member('code-proof.json')
        derived = member('derivation.json')
        transport_join = member('transport-link-proof.json')
        endpoint_join = member('endpoint-link-proof.json')
        for phase in ('finish-build', 'finish-link', 'run'):
            owner = member(phase + '/owner.json')
            assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
            assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
            assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
            assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert finite['finite_component_passed'] and not finite['provider_adopted']
    assert finite['numerical_runs_per_fixture'] == 1
    assert code['passed'] and len(code['objects']) == 3 and code['no_optional_ISA_in_adapter']
    transport = finite['results']['transport']
    for key, expected in dict(case_count=38, active_rows=38, capacity_rows=40,
            packets_per_implementation=5, raw_calls=10, compared_output_words=187960).items():
        assert transport[key] == expected, key
    assert transport['input_unchanged'] and transport['inactive_sentinels_preserved']
    assert transport['full_word_match'] and transport['first_mismatch'] is None
    assert finite['results']['endpoint']['passed']
    assert finite['results']['endpoint']['raw_packets'] == 24
    assert transport_join['passed'] and endpoint_join['passed']
    assert all(row['all_retained_function_opcodes_match'] for proof in (transport_join, endpoint_join)
               for row in proof['objects'])
    assert len(transport_join['declared_GOTPCREL_fields_resolve_exact_targets']) == 2
    provider = root / 'out/native-cpu-v4-raw-provider/provider'
    contract = doc(provider / 'source-contract.json')
    current_path = provider / 'linked-native-proof.json'
    current = doc(current_path)
    objects, parts = contract['native_objects'], contract['adapter_input_objects']
    assert len(objects) == 8 and len(parts) == 7
    assert current['passed'] and [row['object'] for row in current['objects']] == objects
    assert all(row['all_retained_function_opcodes_match'] for row in current['objects'])
    raw = [dict(path=row['candidate']['path'], bytes=row['candidate']['bytes'],
                sha256=row['candidate']['sha256']) for row in code['objects']]
    adapter = code['baseline_adapter']
    math = [row for row in derived['reused_objects'] if Path(row['path']).name in ('dispatch.o', 'v3.o', 'optional.o')]
    expected = {row['path']: {key: row[key] for key in ('bytes', 'sha256')}
                for row in raw + [adapter] + math}
    assert len(expected) == 7
    assert {row['path']: {key: row[key] for key in ('bytes', 'sha256')} for row in parts} == expected
    assert [row['object'] for row in current['raw_v4_objects']] == raw
    for row in current['raw_v4_objects']:
        assert row['all_retained_function_opcodes_match'] and row['functions']
    for name, identity in derived['source_inputs'].items():
        assert pin(root / name) == identity, name
    for row in parts + objects:
        assert pin(root / row['path']) == {key: row[key] for key in ('bytes', 'sha256')}
    inputs = {str(path.relative_to(root)): pin(path) for path in
        (manifest_path, archive, current_path, provider / 'source-contract.json')}
    inputs.update(expected)
    return dict(inputs=inputs, scope=dict(component_result=pin(manifest_path),
        component_link_proof=manifest['members']['out/native-cpu-v4-raw-packets/transport-link-proof.json'],
        component_compile_result=manifest['members']['out/native-cpu-v4-raw-packets/build-result.json'],
        current_AQ_link_proof=pin(current_path), candidate_transport_object=raw[0],
        compared_output_words=187960,
        scope='Sealed original genuine R38/cap40/all187960-word Transport execution and focused Endpoint check, exact same new raw3/baseline adapter and unchanged AO math3. All seven exact parts are joined into current AQ eighth object; all eight provider objects and new raw ranges are joined into the exact library. No fixture rerun. Fresh unchanged C8 must prove actual packet-level v4 selection and DopriPhase execution; no benefit/frame/scientific/release transfer.'))
