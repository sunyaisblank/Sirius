"""Seal the first original C8 gate for the closed three-stage provider."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R = Path.cwd()
A = R / 'attestations/qualification/363aa0c2/native-cpu-closed-all-original-c8'
FOLDERS = ['native-cpu-closed-all-component', 'native-cpu-closed-all-provider', 'native-cpu-closed-all-cohort']


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not A.exists()
owners, sources, products = {}, [Path(__file__).resolve()], []
for folder in FOLDERS:
    for path in (R / 'out' / folder).rglob('*'):
        if not path.is_file() or '__pycache__' in path.parts:
            continue
        (products if path.suffix in ('.so', '.o', '.bc') else sources).append(path)
        if path.name == 'owner.json':
            row = doc(path)
            assert row['status'] == 'terminal' and row['passed'] and row['child_exit'] == 0
            assert row['source_unchanged'] and row['whole_inputs_unchanged']
            assert not row['remaining_owned_processes'] and not row['cleanup_errors']
            assert row['observed_births'] and all(item['absent'] for item in row['observed_births'])
            assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()
            owners[str(path.relative_to(R))] = dict(identity=pin(path), source=row['source'],
                wall_seconds=row['wall_seconds'], peak_sampled_rss_kib=row['peak_sampled_rss_kib'])
assert len(owners) == 4
C = R / 'out/native-cpu-closed-all-component'
P = R / 'out/native-cpu-closed-all-provider/provider'
H = R / 'out/native-cpu-closed-all-cohort'
component, proof, acceptance = doc(C / 'build-result.json'), doc(P / 'linked-native-proof.json'), doc(H / 'control-8-acceptance.json')
assert component['passed'] and component['compiled_interface_passed'] and not component['numerical_qualification_transferred']
assert len(component['parts']) == 11 and len(component['core_definitions']) == 9
assert proof['passed'] and proof['shared_value_interface']['raw_core_calls'] == dict(transport=2696, endpoint=3954, dopri_phase=1576)
assert set(proof['raw_unwind_tables']) == {'transport', 'endpoint', 'dopri_phase'}
for stage, table in proof['raw_unwind_tables'].items():
    assert all(table[key] is True for key in ('full_function_FDE', 'complete_GP_save_rules',
        'shrinkwrap_and_epilogues_checked', 'asynchronous', 'exact_relative_CFI_matches_object'))
    if stage != 'dopri_phase':
        assert table['actual_GP_push_and_restore_machine_bytes_checked']
assert acceptance['passed'] and acceptance['all32_readbacks_match_preserved_originals']
assert acceptance['native_first_packets_handled'] == 3 and acceptance['first_raw_SIMD_entries_completed'] == 2
assert acceptance['observed_Transport12_implementation'] == 'F scalar one-row tasks'
optional = acceptance['optional_math']
assert optional['all_selected_optional_math_observed'] and optional['genuine_private_FP_scope_observed']
assert optional['selected_simd_stages'] == optional['observed_optional_math_stages'] == ['dopri_phase', 'endpoint']
assert optional['component_numerical_qualification_transferred'] is False
products += [R / component['adapter']['path'], R / component['output']['path'], R / proof['library']]
assert pin(R / proof['library']) == proof['library_identity']
predecessor = 'attestations/qualification/dc584b05/native-cpu-closed-remaining-context-native-gates/MANIFEST.json'
assert pin(R / predecessor)['sha256'] == '01e6c6626f6c027b8bfb1d8d053d987334cbefa54d0d37fa21f6e145caaa6c6d'
A.mkdir(parents=True)
bundles = {}
for label, paths in [('source-and-evidence', sources), ('products', products)]:
    entries = {str(path.relative_to(R)):pin(path) for path in sorted(set(paths))}
    target = A / (label + '.zip')
    with zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name in entries:
            archive.write(R / name, name)
    with zipfile.ZipFile(target) as archive:
        assert len(archive.namelist()) == len(entries) and set(archive.namelist()) == set(entries)
        for name, expected in entries.items():
            with archive.open(name) as stream:
                assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected['sha256']
            assert archive.getinfo(name).file_size == expected['bytes'] and pin(R / name) == expected
    assert target.stat().st_size < 100 * 1024 * 1024
    bundles[label] = dict(path=str(target.relative_to(R)), **pin(target), members=entries)
manifest = dict(status='CLOSED_THREE_STAGE_PROVIDER_ORIGINAL_C8_ALL32_PASS_BROADER_GATES_OPEN',
    candidate_revision='363aa0c24c8dc8ba54a1e1e2a22bc4f57675adef',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    owners=owners, predecessor={predecessor:pin(R / predecessor)}, bundles=bundles,
    library=dict(path=proof['library'], **proof['library_identity']),
    original_C8_acceptance=pin(H / 'control-8-acceptance.json'), optional_math=optional,
    C8_observed_scalar_stages=['transport'], C8_observed_closed_SIMD_stages=['endpoint', 'dopri_phase'],
    closed_SIMD_stages_unobserved_by_first_packet_C8=['transport'],
    scope='Actual eleven-part component and nine reused strict local math cores, final exact code/relocation/ISA/FP/live-job/callback joins and full GP/FDE object-to-library proofs for all three closed entries. One original C8 completes all32 full-buffer digest assertions and all3 first native row/pure returns. Its original Transport12 policy remains scalar; only first selected Endpoint/DP SIMD and local-core entry/return paths are passively observed, so no closed Transport SIMD runtime coverage is invented. Counts,51markers54modulejoins, original executable/inputs/assertions/cache/governors and limits are unchanged. C8 debugger timing is not performance evidence. No ABBA22/Decimal/all8, original frame, physical hardware, full scientific/release/main acceptance or qualification transfer; all eight ultimate-goal subgoals remain open.')
with (A / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)), passed_owners=4, complete_readbacks=32,
    observed_closed_SIMD_stages=['endpoint','dopri_phase'], bundles={name:row['bytes'] for name,row in bundles.items()})))
