"""Preserve the metadata-only repair and its original C8 gate, without aliases."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R = Path.cwd()
A = R / 'attestations/qualification/620a0c6e/native-cpu-closed-dp-unwind-original-c8'
FOLDERS = ['native-cpu-closed-dp-unwind', 'native-cpu-closed-dp-unwind-component',
           'native-cpu-closed-dp-unwind-provider', 'native-cpu-closed-dp-unwind-cohort']


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
            assert row['observed_births'] and all(birth['absent'] for birth in row['observed_births'])
            assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()
            owners[str(path.relative_to(R))] = dict(identity=pin(path), source=row['source'],
                wall_seconds=row['wall_seconds'], peak_sampled_rss_kib=row['peak_sampled_rss_kib'])
assert len(owners) == 5
C8 = R / 'out/native-cpu-closed-dp-unwind-cohort/control-8-acceptance.json'
acceptance = doc(C8)
assert acceptance['passed'] and acceptance['all32_readbacks_match_preserved_originals']
assert acceptance['optional_math']['all_selected_optional_math_observed']
proof = doc(R / 'out/native-cpu-closed-dp-unwind-provider/provider/linked-native-proof.json')
assert proof['passed']
library = R / proof['library']
assert pin(library) == proof['library_identity']
products.append(library)
for folder in ('native-cpu-closed-dp-unwind', 'native-cpu-closed-dp-unwind-component'):
    base = R / 'bin/linux-gcc' / folder
    for path in base.rglob('*'):
        if path.is_file():
            (products if path.suffix in ('.bc', '.o', '.so') else sources).append(path)
code = doc(R / 'out/native-cpu-closed-dp-unwind/code-proof.json')
assert code['passed'] and code['full_function_bytes_identical'] and code['normalized_code_relocations_identical']
assert code['unwind_table']['full_function_FDE'] and code['unwind_table']['complete_GP_save_rules']
predecessor = 'attestations/qualification/4332fb95/native-cpu-closed-dp-original-c8-attempt/MANIFEST.json'
assert pin(R / predecessor)['sha256'] == 'f404d0166d95198936ba3d1e32693e303d39db4e10f701d642256ba58ee9bf47'
A.mkdir(parents=True)
bundles = {}
for label, paths in [('source-and-evidence', sources), ('products', products)]:
    entries = {str(path.relative_to(R)): pin(path) for path in sorted(set(paths))}
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
manifest = dict(status='CLOSED_DP_UNWIND_AND_ORIGINAL_C8_PASS_FULL_GOAL_OPEN',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    owners=owners, predecessor={predecessor:pin(R / predecessor)}, bundles=bundles,
    original_C8_acceptance=pin(C8), optional_math=acceptance['optional_math'],
    exact_code_proof=pin(R / 'out/native-cpu-closed-dp-unwind/code-proof.json'),
    scope='The sole private DP uwtable attribute adds asynchronous unwind information while preserving every original instruction byte, original allocated payload and normalized non-unwind relocation. Full-range FDE/CFA/GP-save rules are checked at all calls/returns and joined to the linked provider. Original C8/all32 completes with same-thread positive-job/private-FP-scope optional entries and returns. Previous failed attempt remains independently sealed. This is finite numerical and observation evidence; no performance inference, ABBA22/Decimal/all8, original frame, hardware/scientific/release/main acceptance is transferred. All eight ultimate-goal subgoals remain open; remaining closed Transport/Endpoint source has not been compiled or admitted.')
with (A / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)), passed_owners=len(owners),
    bundles={key:value['bytes'] for key,value in bundles.items()})))
