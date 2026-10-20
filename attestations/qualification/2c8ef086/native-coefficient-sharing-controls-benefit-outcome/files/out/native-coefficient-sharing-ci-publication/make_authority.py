"""Compose actual static eligibility/array authority; never launch its products.

ROOT-only after successful actual publication and PE audit plus a separate ROOT implementation
review. No historical result/hash or runtime right transfers.
"""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
BASE = '2c8ef08693d46f57d5dea14efdbddd4874dac684'
INPUTS = ROOT / 'out/native-coefficient-sharing-producer-inputs'
CHANGED = {'kTransportShader'}


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def pin(path):
    data = path.read_bytes()
    return {'path': str(path.relative_to(ROOT)), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest()}


def value_pin(row):
    return {key: row[key] for key in ('bytes', 'sha256')}


def text(path):
    data = path.read_bytes()
    return data.decode('utf-16' if data.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig')


def local(name):
    path = ROOT / name
    assert not Path(name).is_absolute() and path.resolve(strict=True) == path
    assert path.is_file() and not path.is_symlink()
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('stage', 'revision', 'run-id', 'pe-audit', 'publication-check',
                 'download', 'root-manual-admission'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    assert re.fullmatch('[0-9a-f]{40}', args.revision)
    assert args.run_id.isdigit() and int(args.run_id) > 0
    stage = ROOT / args.stage
    assert args.stage == 'bin/windows-msvc/native-coefficient-sharing-matched-current'
    assert stage.resolve(strict=True) == stage and stage.is_dir()
    eligibility_path = WORK / 'root-actual-producer-eligibility.json'
    authority_path = WORK / 'matched-array-authority.json'
    assert not eligibility_path.exists() and not authority_path.exists()
    paths = {key: local(name) for key, name in (
        ('static_PE_audit', args.pe_audit), ('publication', args.publication_check),
        ('download', args.download), ('root_manual_admission', args.root_manual_admission))}
    audit, publication, download, manual = (read(paths[key]) for key in (
        'static_PE_audit', 'publication', 'download', 'root_manual_admission'))
    assert publication['schema'] == 'sirius-native-coefficient-sharing-actual-publication-check-v1'
    assert publication['passed'] is True and publication['runtime_admitted'] is False
    assert publication['all_source_bytes_match_immutable_Git'] is True
    assert publication['canonical2c8_files_conserved'] == 538
    assert publication['actual_completed_one_attempt_CI_exact'] is True
    assert download['API_digest_CRC_safe_members_and_manifest_exact'] is True
    assert audit['audit_completed'] is True and audit['array_and_private_compile_checks_passed'] is True
    assert audit['full_producer_eligibility_passed'] is False and audit['runtime_admitted'] is False
    assert audit['source_revision'] == publication['producer_source_revision'] == download['producer_source_revision'] == args.revision
    assert str(audit['run_id']) == str(download['run_id']) == args.run_id
    assert set(audit['changed_main_arrays']) == CHANGED
    assert manual['schema'] == 'sirius-native-coefficient-sharing-actual-producer-manual-admission-v1'
    assert manual['passed'] is True and manual['producer_source_revision'] == args.revision
    assert manual['canonical_CXX_source_revision'] == BASE and manual['stage'] == args.stage
    assert manual['run_id'] == args.run_id and manual['runtime_admitted'] is False
    for key in ('static_PE_audit', 'publication', 'download'):
        assert manual[key] == pin(paths[key])
    for key in ('actual_CL5_LINK2_canonical_usage_ROOT_read',
                'canonical_Vulkan_MAP_archive_LINK_selection_ROOT_read',
                'actual_source_tools_SDK_CRT_import_provenance_ROOT_read',
                'original_observer_integrity_reuse_ROOT_validated'):
        assert manual[key] is True
    review_row = manual['separate_ROOT_actual_producer_review']
    review_path = local(review_row['path'])
    assert pin(review_path) == review_row
    review = read(review_path)
    assert review['passed'] is True
    assert review.get('producer_source_revision', review.get('source_revision')) == args.revision

    producer = read(stage / 'producer-result.json')
    manifest = read(stage / 'manifest.json')
    assert producer['producer_profile'] == 'native-coefficient-sharing-one-shader-only'
    assert producer['canonical_CXX_source_revision'] == BASE
    assert type(producer['private_compile_object_count']) is int and producer['private_compile_object_count'] == 5
    assert type(producer['expected_private_compile_objects']) is int and producer['expected_private_compile_objects'] == 5
    assert producer['source_revision'] == producer['github_sha'] == manifest['source_revision'] == args.revision
    assert producer['run_id'] == args.run_id and producer['run_attempt'] == '1'
    assert producer['repository'] == 'sunyaisblank/Sirius' and producer['failure'] is None
    assert all(producer[key] is True for key in ('passed', 'configure_passed', 'build_passed',
                                               'source_unchanged', 'frozen_build_inputs_unchanged'))
    assert all(producer[key] is False for key in ('pe_or_provider_or_test_executed',
        'mandatory_gate_executed', 'native_build_attestation_claimed',
        'scientific_qualification_claimed', 'release_or_adoption_claimed'))
    assert manifest['passed'] is True and manifest['complete_host_or_sdk_or_tool_dso_backup'] is False
    assert manifest['schema'] == 'sirius-native-packed-ci-export-manifest-v1'
    assert producer['schema'] == 'sirius-native-packed-matched-ci-producer-v1'
    manifest_identity = value_pin(pin(stage / 'manifest.json'))
    producer_identity = value_pin(pin(stage / 'producer-result.json'))
    assert download['stage'] == args.stage
    assert manifest_identity == value_pin(download['members']['manifest.json'])
    assert producer_identity == value_pin(download['members']['producer-result.json'])
    assert producer_identity == value_pin(manifest['files']['producer-result.json'])
    before = read(stage / 'build-inputs-before.json')
    assert before == read(stage / 'build-inputs-after.json')
    assert len(before) == publication['before_after_build_rows']
    assert len(read(stage / 'source-before-configure.json')) == publication['all_Git_source_files']
    source_rows = read(stage / 'source-before-configure.json')
    normalize = lambda value: str(value).replace('\\', '/').lower()
    source_root = normalize(producer['source_root']).rstrip('/') + '/'
    relative_source = {normalize(name)[len(source_root):]: row for name, row in source_rows.items()
                       if normalize(name).startswith(source_root)}
    assert len(relative_source) == len(source_rows)
    relative_before = {normalize(name): row for name, row in before.items()}

    assert len(publication['selected_tools']) == 15
    assert publication['PE_imports']['baseline'] == publication['PE_imports']['candidate']
    packaging = read(stage / 'runtime-packaging.json')
    assert publication['runtime_packaging_versions'] == [packaging['msvc_version'], packaging['release_crt_version']]
    compiler_minor = re.fullmatch(r'19\.(\d+)\..+', packaging['msvc_version'])
    assert compiler_minor and packaging['release_crt_version'].startswith('14.' + compiler_minor[1] + '.')

    observer = INPUTS / 'observer/audited-observer.cpp'
    assert value_pin(pin(observer)) == {'bytes': 342470,
        'sha256': 'c9485d03ec55b5c1307a34097f15efbc58c4f84749e9f3736ad063a68676b601'}
    test = 'tests/backend/retained_compute_test.cpp'
    old_test_blob = subprocess.check_output(['git', 'rev-parse',
        '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e:' + test], cwd=ROOT).decode().strip()
    assert subprocess.check_output(['git', 'rev-parse', BASE + ':' + test], cwd=ROOT).decode().strip() == old_test_blob
    composition = read(INPUTS / 'root-composed-inputs.json')
    assert composition['schema'] == 'sirius-native-coefficient-sharing-diagnostic-input-composition-v1'
    assert composition['passed'] is True and composition['source_revision'] == BASE
    assert len(composition['files']) == 77 and composition['excluded_SPVs_exact'] == 34
    assert composition['main_arrays'] == {'changed': 1, 'exact': 39, 'seven_tables_and_other_text_exact': True}
    for name, row in composition['files'].items():
        path = INPUTS / name
        assert path.resolve(strict=True) == path and path.is_file() and not path.is_symlink()
        assert value_pin(pin(path)) == value_pin(row)
        suffix = 'out/native-coefficient-sharing-producer-inputs/' + name
        assert value_pin(relative_source[suffix.lower()]) == value_pin(row)
        assert value_pin(relative_before[source_root + suffix.lower()]) == value_pin(row)

    assert value_pin(pin(INPUTS / 'root-composed-inputs.json')) == value_pin(relative_source['out/native-coefficient-sharing-producer-inputs/root-composed-inputs.json'])
    assert len(audit['private_compiles']) == 5
    for record in audit['private_compiles'].values():
        assert len(record['command_tlogs']) >= 1 and len(record['read_tlogs']) >= 1
        command = '\n'.join(text(stage / name) for name in record['command_tlogs']).lower()
        for flag in ('/w4', '/wx', '/o2', '/ob2', '/md', '/fp:precise', '/std:c++latest',
                     '/permissive-', 'sirius_contract_mode=2'):
            assert flag in command
    canonical = audit['canonical_device']
    assert canonical['actual_CL_source_records_joined'] is True
    backend = stage / 'metadata/canonical-backend/src/sirius/backend'
    assert value_pin(pin(backend / 'Release/sirius_backend.lib')) == canonical['archive']
    device_paths = [path for path in (backend / 'sirius_backend.dir/Release').rglob('*.obj')
                    if path.name.lower() in ('vulkan_device.obj', 'vulkan_device.cpp.obj')]
    assert len(device_paths) == 1 and value_pin(pin(device_paths[0])) == canonical['object']
    assert canonical['exact_archive_member_label']

    arms = {}
    for name in ('baseline', 'candidate'):
        product = audit['products'][name]
        assert len(product['main_arrays']) == 40 and len(product['owner_arrays']) == 2
        exe = stage / 'payload' / name / ('original_' + name + '_controls.exe')
        assert value_pin(pin(exe)) == product['PE']
        raw = exe.read_bytes()
        mapping_path = exe.with_suffix('.map')
        assert value_pin(pin(mapping_path)) == product['map']
        mapping = text(mapping_path)
        device_members = re.findall(r'sirius_backend(?:\.lib)?[:(][^\s)]*vulkan_device\.obj', mapping, re.I)
        assert len({value.lower() for value in device_members}) == 1
        assert sorted(set(device_members)) == product['canonical_device_map_members']
        assert not re.search(r'sirius_backend(?:\.lib)?[:(][^\s)]*retained_compute\.obj', mapping, re.I)
        for category in ('main_arrays', 'owner_arrays'):
            for row in product[category].values():
                assert row['locations']
                for location in row['locations']:
                    offset = location['file_offset']; count = row['bytes']
                    assert 0 <= offset <= len(raw) - count and count > 0
                    assert hashlib.sha256(raw[offset:offset + count]).hexdigest() == row['sha256']
        available = {}
        header = INPUTS / name / 'retained_kernels.h'
        for role in ('transport', 'endpoint'):
            values = re.findall(r'inline constexpr bool k' + role.title() + r'FmaAvailable = (true|false);', header.read_text())
            assert len(values) == 1; available[role] = values[0] == 'true'
        arms[name] = {'PE': product['PE'], 'arrays': product['main_arrays'],
                      'owner_arrays': product['owner_arrays'], 'fma_available': available}
        links = product['link_command_tlogs']; assert len(links) == 1
        command = '\n'.join(line for line in text(stage / links[0]).splitlines() if not line.startswith('^'))
        assert '/verbose:lib' in command.lower()
        assert not re.search(r'(?i)(?<!\S)/(?:force|wholearchive)(?=[:\s]|$)', command)
        reads = product['link_read_tlogs']; assert len(reads) == 1
        read_inputs = {normalize(line.strip()) for line in text(stage / reads[0]).splitlines()
                       if line.strip() and not line.startswith('^')}
        canonical_archive = normalize(producer['build_root']).rstrip('/') + '/src/sirius/backend/release/sirius_backend.lib'
        assert canonical_archive in read_inputs
        normalized = re.sub('baseline|candidate', '<ARM>', command, flags=re.I)
        if name == 'baseline': baseline_link = normalized
        else: assert normalized == baseline_link
        for dll, origin in packaging['runtime_origins'].items():
            assert value_pin(pin(stage / 'payload' / name / dll)) == value_pin(origin)
        imports_path = paths['publication'].parent / (name + '-PE-imports.stdout')
        assert re.findall(r'DLL Name:\s*(\S+)', imports_path.read_text()) == publication['PE_imports'][name]
    left, right = arms['baseline']['arrays'], arms['candidate']['arrays']
    assert set(left) == set(right) and {key for key in left if value_pin(left[key]) != value_pin(right[key])} == CHANGED
    tables = {key for key in left if key.endswith('Program')}
    assert len(tables) == 7 and all(value_pin(left[key]) == value_pin(right[key]) for key in tables)
    assert {key: value_pin(row) for key, row in arms['baseline']['owner_arrays'].items()} == {key: value_pin(row) for key, row in arms['candidate']['owner_arrays'].items()}
    assert arms['baseline']['fma_available'] == arms['candidate']['fma_available']
    for header in ('retained_owner_kernels.h', 'retained_owner_layout.h'):
        assert (INPUTS / 'baseline' / header).read_bytes() == (INPUTS / 'candidate' / header).read_bytes()

    scope = ('One successful matched Windows compile-only producer: five private strict objects/two links, '
             'one changed native32 array and34 excludedSPVs, all40main+2owner payloads/7tables, canonical '
             'Vulkan archive selection, actual bounded source/tool/fetched/configured/SDK+CRT/import joins. '
             'ROOT manual and separate ROOT actual-producer reviews remain bound evidence, not a full-host '
             'backup or numerical/benefit/science/Mandatory/release/adoption/runtime right.')
    eligibility = {'schema': 'sirius-native-coefficient-sharing-actual-producer-eligibility-v1', 'passed': True,
        'producer_source_revision': args.revision, 'canonical_CXX_source_revision': BASE,
        'run_id': args.run_id, 'stage': args.stage, **{key: pin(path) for key, path in paths.items()},
        'separate_ROOT_actual_producer_review': pin(review_path), 'five_private_actual_strict_compile_records_joined': True,
        'two_LINK_commands_equal_modulo_arm_names': True, 'canonical_device': canonical,
        'manifest': pin(stage / 'manifest.json'), 'producer': pin(stage / 'producer-result.json'),
        'canonical_Vulkan_selection_evidence': 'Exact unique archive object, four Vulkan function MAP member joins, actual LINK.read and LINK command. /VERBOSE:LIB requested; saved stdout has no individual Loaded extraction lines, so none are asserted.',
        'runtime_packaging': pin(stage / 'runtime-packaging.json'),
        'original_observer_integrity_reuse': {'complete_observer': pin(observer),
            'unchanged_canonical_test_blob': old_test_blob,
            'prior_complete_matching_integrity': 'fd530d957d2e0737f5453d88d78d979f8ed17607'},
        'warning_scope': publication['warning_codes_including_summaries'],
        'runtime_admitted': False, 'scope': scope}
    eligibility_bytes = (json.dumps(eligibility, indent=2) + '\n').encode()
    authority = {'schema': 'sirius-native-packed-matched-array-authority-v1', 'passed': True,
        'producer_source_revision': args.revision, 'canonical_CXX_source_revision': BASE,
        'actual_PE_payload_and_selector_joins_passed': True,
        'original27_test_bodies_and1492_assertion_lines_equal': True,
        'seven_program_tables_and_public_constants_equal': True, 'owner2_payloads_equal': True,
        'manifest_identity': manifest_identity, 'producer_identity': producer_identity,
        'arms': arms, 'root_actual_producer_eligibility': {'path': str(eligibility_path.relative_to(ROOT)),
            'bytes': len(eligibility_bytes), 'sha256': hashlib.sha256(eligibility_bytes).hexdigest()},
        'runtime_admitted': False, 'scope': scope}
    authority_bytes = (json.dumps(authority, indent=2) + '\n').encode()
    # All data checks precede exclusive output writes; no output enters its own inputs.
    for path, data in ((eligibility_path, eligibility_bytes), (authority_path, authority_bytes)):
        with path.open('xb') as stream: stream.write(data)
    print(json.dumps({'producer_eligible': True, 'matched_array_authority_written': True,
                      'runtime_admitted': False}))


if __name__ == '__main__':
    main()
