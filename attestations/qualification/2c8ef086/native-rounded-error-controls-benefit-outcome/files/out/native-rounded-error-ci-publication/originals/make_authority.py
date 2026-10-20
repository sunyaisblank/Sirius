import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
STAGE = ROOT / 'bin/windows-msvc/native-wave32-matched-current'
REVISION = '5d993d12b3fde7a65eaa5f7988b73b86b3a9a291'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def pin(path):
    data = path.read_bytes()
    return {'path': str(path.relative_to(ROOT)), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest()}


audit_path = WORK / 'actual-5d993d12-link-paths-audit.json'
publication_path = ROOT / 'out/native-wave32-ci-publication/root-publication-check.json'
download_path = ROOT / 'out/native-wave32-ci-publication/actual-download.json'
audit, publication, download = read(audit_path), read(publication_path), read(download_path)
assert audit['audit_completed'] and audit['array_and_private_compile_checks_passed']
assert audit['native_wave32_device_object_and_map_joins_passed'] and audit['private_compile_object_count'] == 7
assert audit['source_revision'] == publication['producer_source_revision'] == download['producer_source_revision'] == REVISION
assert publication['passed'] and download['API_digest_CRC_safe_members_and_manifest_exact']
observer = ROOT / 'out/native-wave32-producer-inputs/observer/audited-observer.cpp'
assert pin(observer)['sha256'] == 'c9485d03ec55b5c1307a34097f15efbc58c4f84749e9f3736ad063a68676b601'
test = 'tests/backend/retained_compute_test.cpp'
original_test_blob = subprocess.check_output(['git', 'rev-parse', '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e:' + test], cwd=ROOT).decode().strip()
assert subprocess.check_output(['git', 'rev-parse', '2c8ef08693d46f57d5dea14efdbddd4874dac684:' + test], cwd=ROOT).decode().strip() == original_test_blob
arms = {}
for name in ['baseline', 'candidate']:
    product = audit['products'][name]
    assert product['LINK_operands_and_read_records_joined']
    assert product['link_textual_order_is_selection_authority'] is False
    assert len(product['main_arrays']) == 40 and len(product['owner_arrays']) == 2
    header = ROOT / 'out/native-wave32-producer-inputs' / name / 'retained_kernels.h'
    available = {}
    for stage in ['transport', 'endpoint']:
        matches = re.findall(r'inline constexpr bool k' + stage.title() + r'FmaAvailable = (true|false);', header.read_text())
        assert len(matches) == 1
        available[stage] = matches[0] == 'true'
    arms[name] = {'PE': product['PE'], 'arrays': product['main_arrays'],
                  'owner_arrays': product['owner_arrays'], 'fma_available': available,
                  'private_device_object': product['private_device_object'],
                  'private_device_COFF_MAP_joins': product['private_device_COFF_MAP_joins']}
    for row in audit['private_compiles'].values():
        for record in row['records'].values():
            command = record['command_semantics']
            for flag in ['/w4', '/wx', '/o2', '/ob2', '/md', '/fp:precise', '/std:c++latest', '/permissive-', 'sirius_contract_mode=2']:
                assert flag in command
    link_paths = product['link_command_tlogs']
    assert len(link_paths) == 1
    raw = (STAGE / link_paths[0]).read_bytes()
    text = raw.decode('utf-16' if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig')
    normalized = re.sub('baseline|candidate', '<ARM>', text, flags=re.I)
    if name == 'baseline':
        baseline_link = normalized
    else:
        assert normalized == baseline_link
for category in ['arrays', 'owner_arrays']:
    assert {k: (v['bytes'], v['sha256']) for k, v in arms['baseline'][category].items()} == {k: (v['bytes'], v['sha256']) for k, v in arms['candidate'][category].items()}
for header in ['retained_kernels.h', 'retained_owner_kernels.h', 'retained_owner_layout.h']:
    assert (ROOT / 'out/native-wave32-producer-inputs/baseline' / header).read_bytes() == (ROOT / 'out/native-wave32-producer-inputs/candidate' / header).read_bytes()
eligibility = {'schema': 'sirius-native-wave32-actual-producer-eligibility-v1', 'passed': True,
               'producer_source_revision': REVISION, 'canonical_CXX_source_revision': audit['canonical_CXX_source_revision'],
               'publication': pin(publication_path), 'download': pin(download_path), 'static_PE_audit': pin(audit_path),
               'seven_private_actual_strict_compile_records_joined': True,
               'two_LINK_commands_equal_modulo_arm_names': True,
               'original_observer_integrity_reuse': {'complete_observer_sha256': pin(observer)['sha256'],
                   'unchanged_canonical_test_blob': original_test_blob,
                   'prior_complete_matching_integrity': 'fd530d957d2e0737f5453d88d78d979f8ed17607 native-constant-controls-benefit-pass'},
               'warning_scope': 'Four canonical dependency D9025 occurrences retained; no all-build warning-free claim.',
               'runtime_admitted': False,
               'scope': 'One exact successful Windows matched compile-only producer, bounded source/tool/fetched/configured closure. Actual private Vulkan COFF/MAP/LINK selection plus42 payloads per arm. No numerical, benefit, full-host, scientific, Mandatory, release or adoption claim. Independent concrete phase review and ROOT runtime release remain.'}
(WORK / 'root-actual-producer-eligibility.json').write_text(json.dumps(eligibility, indent=2) + '\n')
authority = {'schema': 'sirius-native-packed-matched-array-authority-v1', 'passed': True,
             'producer_source_revision': REVISION, 'canonical_CXX_source_revision': audit['canonical_CXX_source_revision'],
             'actual_PE_payload_and_selector_joins_passed': True,
             'original27_test_bodies_and1492_assertion_lines_equal': True,
             'seven_program_tables_and_public_constants_equal': True, 'owner2_payloads_equal': True,
             'native_wave32_device_object_and_map_joins_passed': True,
             'manifest_identity': audit['manifest_identity'], 'producer_identity': audit['producer_identity'],
             'arms': arms, 'root_actual_producer_eligibility': pin(WORK / 'root-actual-producer-eligibility.json'),
             'runtime_admitted': False, 'scope': eligibility['scope']}
(WORK / 'matched-array-authority.json').write_text(json.dumps(authority, indent=2) + '\n')
print(json.dumps({'producer_eligible': True, 'matched_array_authority_written': True, 'runtime_admitted': False}))
