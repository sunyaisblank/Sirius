"""Read sealed packed/raw and arithmetic evidence without restoring retired products."""
from pathlib import Path
import hashlib
import json
import re
import zipfile

R = Path.cwd()
MANIFEST = R / 'attestations/qualification/e17097b9/native-cpu-inline-private-bank-finite-compatibility/MANIFEST.json'
MANIFEST_SHA = '757366c154706de8085e8c5c00be7d0fc99dfca8faefd8905905c5d2b7ec74a8'
MATH = R / 'attestations/qualification/ca334e6e/native-cpu-exact-product-zero-residual-finite-compatibility/MANIFEST.json'
MATH_SHA = '44d228597b02d724ad7a5b0ebfe8d69073bf306977bacd3a93013d1462213a66'
AQ = R / 'out/native-cpu-packed-math-provider/provider'
STAGES = ('transport', 'endpoint', 'dopri_phase')


def doc(path):
    return json.loads(path.read_text())


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def checked(path, row):
    assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}, path


def terminal(path):
    value = doc(path)
    assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    return value


def contract_inputs(contract):
    for name, expected in contract['source_inputs'].items():
        checked(R / name, expected)


def bundle(path, sha):
    assert isinstance(sha, str) and re.fullmatch(r'[0-9a-f]{64}', sha), 'AS component manifest is not yet bound'
    assert pin(path)['sha256'] == sha
    manifest = doc(path)
    assert manifest['finite_component_passed'] and not manifest['provider_adopted']
    archive = R / manifest['archive']['path']
    checked(archive, manifest['archive'])
    with zipfile.ZipFile(archive) as zipped:
        assert set(zipped.namelist()) == set(manifest['members'])
        files = {name: zipped.read(name) for name in zipped.namelist()}
    for name, expected in manifest['members'].items():
        assert data_pin(files[name]) == expected, name
    for row in manifest['owners']:
        checked_bytes = files[row['path']]
        assert data_pin(checked_bytes) == {key: row[key] for key in ('bytes', 'sha256')}
        value = json.loads(checked_bytes)
        assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
        assert value['source'] == row['source'] and value['source_unchanged'] and value['whole_inputs_unchanged']
        assert not value['remaining_owned_processes'] and not value['cleanup_errors']
        assert value['observed_births'] and all(item['absent'] for item in value['observed_births'])
        assert row['all_owned_births_reaped']
        folder = row['path'].removesuffix('owner.json')
        assert files[folder + 'inputs-before.json'] == files[folder + 'inputs-after.json']
    return manifest, archive, files


def component():
    manifest, archive, files = bundle(MANIFEST, MANIFEST_SHA)
    math, math_archive, math_files = bundle(MATH, MATH_SHA)
    assert manifest['math_component_manifest'] == dict(path=str(MATH.relative_to(R)), **pin(MATH))
    assert manifest['new_raw_compilations'] == 3 and manifest['new_CXX_compilations'] == 0
    assert manifest['fixture_runs_per_fixture'] == 1
    assert manifest['original_v3_dispatcher_AQ_baseline_adapter_and_first7_AP_native_objects_reused']
    assert manifest['original_fixture_objects_restored_exactly_from_sealed_archive']
    assert manifest['all9380_literal_operators_and_original_masks_preserved']
    assert manifest['only3_opcode_inline_qualifiers_and48_integer_privateD_zero_words']
    prefix = 'out/native-cpu-inline-private-bank/'
    read = lambda name: json.loads(files[prefix + name])
    derivation, build, linked, compatibility, code = [read(name) for name in (
        'derivation.json', 'build-result.json', 'link-result.json', 'compatibility.json', 'code-proof.json')]
    assert build['passed'] and build['new_raw_compilations'] == 3 and build['new_CXX_compilations'] == 0
    assert linked['passed'] and linked['new_CXX_compilations'] == 0 and compatibility['passed']
    assert code == manifest['code'] and code['passed'] and len(code['objects']) == 3
    assert code['explicit_compiler_private_temporary_bytes_per_store'] == 576
    assert code['removed_EvaluateOneStatic_functions'] == 16
    assert code['all_executable_sections_checked'] and code['unknown_or_non_v4_forms_rejected']
    assert code['new_ISA_vocabulary_admitted'] == []
    assert derivation['literal_instruction_calls_preserved'] == 9380
    assert derivation['math_component_manifest'] == dict(path=str(MATH.relative_to(R)), **pin(MATH))
    assert derivation['compiler_private_storage'] == manifest['compiler_private_storage']
    assert derivation['parent_component_manifest'] == manifest['parent_component_manifest']
    assert derivation['compiler_private_storage'] == dict(payload_words=120, tag_words=24, bytes=576,
        slots=dict(destination=0, a=1, b=2), physical_lanes=8, stack_reuse_established=False,
        host_allocations_changed=False, KernelContext_layout_changed=False)
    aq, proof = doc(AQ / 'source-contract.json'), doc(AQ / 'linked-native-proof.json')
    assert proof['passed'] and derivation['reused_objects'][:7] == aq['native_objects'][:7]
    raw = []
    for slug, row, product in zip(STAGES, code['objects'], build['objects']):
        assert row['stage'] == slug and row['candidate'] == product
        assert row['original'] == proof['raw_v4_objects'][len(raw)]['object']
        assert row['no_FMA'] and row['no_x87'] and row['no_FP_control_writes']
        assert row['EvaluateOneStatic_symbols_calls_absent']
        assert row['all_executable_bytes_decoded'] and row['only_sealed_AR_instruction_forms']
        assert row['new_mnemonics'] == row['new_encoding_forms'] == []
        assert all(name in {'memset', 'sirius_simd_fmaxf'} or name.startswith('sirius_native_v4_')
                   for name in row['undefined'])
        checked(R / product['path'], product)
        assert {key: product[key] for key in ('bytes', 'sha256')} == manifest['products'][product['path']]
        stage = derivation['stages'][slug]
        assert stage['whole_source_inverse_restores_parent']
        assert stage['sole_opcode_inline_qualifier'] and stage['helper_include_rebound']
        assert stage['layer_scopes_unchanged'] == {'transport':159,'endpoint':195,'dopri_phase':45}[slug]
        data = files[stage['derived_path']]
        assert data_pin(data) == derivation['derived'][stage['derived_path']]
        calls = re.findall(rb'^\s+(sirius_native_(?:dense|tagged)_\w+8)\(([^;]*?)\);', data, re.M)
        assert len(calls) == 9
        for name, args in calls:
            assert args.startswith(b'SiriusBank, 2, 0, 1,') if b'dense' in name else args.startswith(b'SiriusPackedBank, SiriusPackedTags,')
        assert not re.search(rb'sirius_native_(?:dense|tagged)_\w+8\([^;]*SiriusContext', data)
        assert manifest['source_math_call_boundaries'][slug] == dict(dense_private_calls=4,
            tagged_private_calls=5, full_original_bank_passed_to_CXX_math=False)
        raw.append(product)
    results = {row['stage']: row['result'] for row in compatibility['fixtures']}
    assert list(results) == ['transport', 'endpoint'] and results == manifest['results']
    for stage in ('transport', 'endpoint'):
        output, association = read(stage + '-output.json'), read(stage + '-link-proof.json')
        assert not output['stderr'] and json.loads(output['stdout']) == results[stage]
        assert association['passed'] and all(row['all_retained_function_opcodes_match'] for row in association['objects'])
        assert manifest['complete_link_association'][stage] == dict(objects=len(association['objects']),
            functions=sum(len(row['functions']) for row in association['objects']),
            declared_GOT_fields=len(association['declared_GOTPCREL_fields_resolve_exact_targets']))
    t, e = results['transport'], results['endpoint']
    for name, expected in dict(passed=True, case_count=38, active_rows=38, capacity_rows=40,
            packets_per_implementation=5, raw_calls=10, compared_output_words=187960,
            mismatch_count=0, first_mismatch=None).items():
        assert t[name] == expected, name
    assert t['input_unchanged'] and t['inactive_sentinels_preserved'] and t['full_word_match']
    assert files[prefix + 'run/transport/input-before.bin'] == files[prefix + 'run/transport/input-after.bin']
    assert files[prefix + 'run/transport/baseline.bin'] == files[prefix + 'run/transport/candidate.bin']
    assert len(files[prefix + 'run/transport/candidate.bin']) == 187960 * 4
    for name, expected in dict(passed=True, selection_rows=23, capacity=32, header_values=6, layouts=2,
            reference_rows=138, scalar_packets=360, raw_packets=24, pure_calls=21199,
            other_program_word_refusals=21142, boundary_refusals=21,
            compared_output_words=5189312960, compared_input_words=601713080).items():
        assert e[name] == expected, name
    for name in ('whole3829_words_and_inactive_capacity_and_canaries', 'input_immutable',
                 'MXCSR_controls_unchanged', 'caller_MXCSR_state_restored'):
        assert e[name] is True, name
    mp = 'out/native-cpu-exact-product-zero-residual/'
    mr = lambda name: json.loads(math_files[mp + name])
    compiled, math_code, compared, primitive = [mr(name) for name in (
        'compile-result.json', 'code-proof.json', 'compatibility.json', 'primitive-compatibility.json')]
    optional = manifest['new_optional_math']
    assert optional in derivation['reused_objects']
    assert optional['path'] == 'bin/linux-gcc/native-cpu-exact-product-zero-residual/optional.o'
    checked(R / optional['path'], optional)
    optional_pin = {key: optional[key] for key in ('bytes', 'sha256')}
    assert optional_pin == math['products'][optional['path']] == compiled['optional_object'] == math_code['optional_object'] == compared['optional_object']
    assert compiled['exit'] == 0 and not compiled['driver_recompiled']
    assert not compiled['original_v3_recompiled'] and not compiled['baseline_dispatcher_recompiled']
    assert math_code['passed'] and math_code['no_FMA'] and math_code['original_v3_and_baseline_dispatcher_byte_preserved']
    assert math_code['wide_product_sites'] == math_code['retention_multiply_sites'] == len(math_code['retention_dependencies']) == 36
    assert math['code']['full_proof'] == dict(member=mp + 'code-proof.json', **data_pin(math_files[mp + 'code-proof.json']))
    assert compared['passed'] and compared['comparison'] == json.loads(math_files[mp + 'run/stdout.log'])
    c = compared['comparison']
    for name, expected in dict(passed=True, calls=555660, compared_words=66679200,
            tag_words=13335840, active_lane_results=1666980, case_tuples=49,
            input_alias_and_inactive_words_checked=True, FP_controls_unchanged=True).items():
        assert c[name] == expected, name
    assert primitive['passed'] and primitive['result'] == json.loads(math_files[mp + 'primitive-run/stdout.log'])
    pr = primitive['result']
    for name, expected in dict(passed=True, case_count=24, candidate_calls=24, historical_primitive_calls=24,
            compared_hi_low_words=384, compared_canary_words=384, independent_center_cases=19,
            denormal_reference_witness=True).items():
        assert pr[name] == expected, name
    assert len(pr['cases']) == 24
    for row in pr['cases']:
        assert row['passed'] and row['whole_words_equal'] and row['center_and_zero_signs_equal'] and row['caller_restored']
        assert row['old_MXCSR'] == row['new_MXCSR']
    oldparts = aq['adapter_input_objects']
    adapter, = [row for row in oldparts if Path(row['path']).name == 'rows_adapter.o']
    preserved = [row for row in oldparts if Path(row['path']).name in ('v3.o', 'dispatch.o')]
    assert len(preserved) == 2
    parts = [adapter] + raw + preserved + [optional]
    for row in parts:
        checked(R / row['path'], row)
    transport = read('transport-link-proof.json')
    joined = {str(Path(row['object']['path']).resolve()): row['object'] for row in transport['objects']}
    for row in parts + aq['native_objects'][:7]:
        assert {key: joined[str((R / row['path']).resolve())][key] for key in ('bytes', 'sha256')} == {key: row[key] for key in ('bytes', 'sha256')}
    return dict(parts=parts, raw_objects=raw, adapter=adapter, manifest=MANIFEST, archive=archive,
        math_manifest=MATH, math_archive=math_archive, math_code=math_code,
        math_code_proof=math['code']['full_proof'], code=code, finite_results=results)
