"""Consume the sealed raw qualification directly; retired phase paths stay archived."""
from pathlib import Path
import hashlib
import json
import zipfile

R = Path.cwd()
MANIFEST = R / 'attestations/qualification/984ba457/native-cpu-v4-raw-packets-finite-compatibility/MANIFEST.json'
MANIFEST_SHA = 'b41bc17d44c4757c51415a1282da22692c85cbd589e92c11724e47da01cdc19b'
AP = R / 'out/native-cpu-endpoint-header-provider/provider'
STAGES = ('transport', 'endpoint', 'dopri_phase')


def doc(path):
    return json.loads(path.read_text())


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def checked(path, expected):
    assert pin(path) == {key: expected[key] for key in ('bytes', 'sha256')}, path


def terminal(path):
    value = doc(path)
    assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    return value


def contract_inputs(contract):
    for name, expected in contract['source_inputs'].items():
        checked(R / name, expected)


def component():
    assert pin(MANIFEST)['sha256'] == MANIFEST_SHA
    manifest = doc(MANIFEST)
    assert manifest['finite_component_passed'] and not manifest['provider_adopted']
    assert manifest['compilation_runs'] == manifest['fixture_runs_per_fixture'] == 1
    assert manifest['new_raw_stage_compilations'] == 3
    assert manifest['new_baseline_adapter_and_fixture_guard_compilations'] == 2
    assert manifest['original_math_raw3_scalar_admission_and_fixture_compilations'] == 0
    archive = R / manifest['archive']['path']
    checked(archive, manifest['archive'])
    with zipfile.ZipFile(archive) as zipped:
        assert set(zipped.namelist()) == set(manifest['members'])
        for name, expected in manifest['members'].items():
            assert data_pin(zipped.read(name)) == expected, name
        def member(name):
            return json.loads(zipped.read('out/native-cpu-v4-raw-packets/' + name))
        derivation, build, code = (member(name) for name in ('derivation.json', 'build-result.json', 'code-proof.json'))
        finite, compatibility, linked = (member(name) for name in ('finite-result.json', 'compatibility.json', 'link-result.json'))
        assert finite['finite_component_passed'] and not finite['provider_adopted']
        assert finite['compilation_runs'] == finite['numerical_runs_per_fixture'] == 1
        assert not finite['original_raw3_math3_and_fixtures_recompiled']
        assert finite['original_failed_setup_phases_preserved'] == ['build', 'link']
        for phase in ('build', 'finish-build', 'link', 'finish-link', 'run'):
            value = member(phase + '/owner.json')
            assert value['status'] == 'terminal'
            assert value['passed'] == (phase not in ('build', 'link'))
            assert value['child_exit'] == (1 if phase in ('build', 'link') else 0)
            assert value['source_unchanged'] and value['whole_inputs_unchanged']
            assert not value['remaining_owned_processes'] and not value['cleanup_errors']
            assert value['observed_births'] and all(row['absent'] for row in value['observed_births'])
            assert zipped.read('out/native-cpu-v4-raw-packets/' + phase + '/inputs-before.json') == zipped.read(
                'out/native-cpu-v4-raw-packets/' + phase + '/inputs-after.json')
            declared, = [row for row in manifest['owners'] if row['path'].endswith('/' + phase + '/owner.json')]
            assert declared['passed'] == value['passed'] and declared['source'] == value['source']
            assert declared['all_owned_births_reaped']
        assert build['passed'] and not build['original_compile_phase_passed']
        assert build['fixture_copy_verification_completed'] and build['new_compilations'] == 0
        assert not build['original_math_and_raw3_recompiled']
        for row in build['products']:
            assert {key: row[key] for key in ('bytes', 'sha256')} == manifest['products'][row['path']]
        assert len(build['fixture_reference_renames']) == 2
        for row in build['fixture_reference_renames']:
            assert row['all_allocated_sections_identical'] and row['all_symbols_and_relocations_inverse_normalized_equal']
        assert compatibility['passed'] and compatibility['fixtures'] == [dict(stage=s, exit=0) for s in STAGES[:2]]
        assert linked['passed'] and finite['results'] == manifest['results']
        for stage in STAGES[:2]:
            output = member(stage + '-output.json')
            assert not output['stderr'] and json.loads(output['stdout']) == finite['results'][stage]
            proof = member(stage + '-link-proof.json')
            assert proof['passed'] and all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
            summary = dict(objects=len(proof['objects']), functions=sum(len(row['functions']) for row in proof['objects']),
                           exact_GOT_fields=len(proof['declared_GOTPCREL_fields_resolve_exact_targets']))
            assert summary == finite['complete_link_association'][stage] == manifest['complete_link_association'][stage]
        t, e = finite['results']['transport'], finite['results']['endpoint']
        for key, expected in dict(passed=True, case_count=38, active_rows=38, capacity_rows=40,
                packets_per_implementation=5, raw_calls=10, compared_output_words=187960,
                mismatch_count=0, first_mismatch=None).items():
            assert t[key] == expected, key
        assert t['input_unchanged'] and t['inactive_sentinels_preserved'] and t['full_word_match']
        prefix = 'out/native-cpu-v4-raw-packets/run/transport/'
        assert zipped.read(prefix + 'input-before.bin') == zipped.read(prefix + 'input-after.bin')
        assert zipped.read(prefix + 'baseline.bin') == zipped.read(prefix + 'candidate.bin')
        assert len(zipped.read(prefix + 'candidate.bin')) == 187960 * 4
        for key, expected in dict(passed=True, selection_rows=23, capacity=32, header_values=6, layouts=2,
                reference_rows=138, scalar_packets=360, raw_packets=24, pure_calls=21199,
                other_program_word_refusals=21142, boundary_refusals=21,
                compared_output_words=5189312960, compared_input_words=601713080).items():
            assert e[key] == expected, key
        for key in ('whole3829_words_and_inactive_capacity_and_canaries', 'input_immutable',
                    'MXCSR_controls_unchanged', 'caller_MXCSR_state_restored'):
            assert e[key] is True, key
        assert code == finite['code'] == manifest['code'] and code['passed']
        assert code['no_optional_ISA_in_adapter'] and code['original_scalar_v3_and_admission_edges_present']
        assert code['all_emitted_instruction_mnemonics_reviewed_against_original_and_explicit_v4_additions']
        assert code['arithmetic_helpers_unchanged'] and len(code['objects']) == 3
        assert derivation['literal_instruction_calls_preserved'] == 9380
        assert derivation['adapter_inverse_restores_original']
        ap = doc(AP / 'source-contract.json')
        assert derivation['reused_objects'][:7] == ap['native_objects'][:7]
        raw = []
        for slug, row in zip(STAGES, code['objects']):
            assert row['stage'] == slug and row['original'] == ap['native_objects'][4 + len(raw)]
            assert row['no_FMA'] and row['direct_optional_helpers_only']
            assert row['candidate']['path'] == 'bin/linux-gcc/native-cpu-v4-raw-packets/objects/' + slug + '.o'
            source = derivation['stages'][slug]
            data = zipped.read(source['derived_path'])
            assert data_pin(data) == derivation['derived'][source['derived_path']]
            assert data.startswith(source['prefix'].encode())
            assert data_pin(data[len(source['prefix'].encode()):]) == source['basis']
            assert source['inverse_restores_exact_AN_source']
            checked(R / row['candidate']['path'], row['candidate'])
            assert {key: row['candidate'][key] for key in ('bytes', 'sha256')} == manifest['products'][row['candidate']['path']]
            raw.append(row['candidate'])
        adapter = code['baseline_adapter']
        checked(R / adapter['path'], adapter)
        assert {key: adapter[key] for key in ('bytes', 'sha256')} == manifest['products'][adapter['path']]
        math = [row for row in ap['adapter_input_objects'] if Path(row['path']).name in ('dispatch.o', 'v3.o', 'optional.o')]
        assert len(math) == 3 and derivation['reused_objects'][7:] == math
        for row in math:
            checked(R / row['path'], row)
    return dict(parts=[adapter] + raw + math, raw_objects=raw, adapter=adapter,
        manifest=MANIFEST, archive=archive, code=code,
        finite_proof=manifest['members']['out/native-cpu-v4-raw-packets/finite-result.json'])
