"""Seal the closed DP compiler, retained-object audits and original C8 evidence."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R = Path.cwd()
V = R / 'out/native-cpu-strict-ir-linking'
A = R / 'attestations/qualification/4332fb95/native-cpu-closed-dp-original-c8-attempt'
FOLDERS = ['native-cpu-canonical-closed-dp', 'native-cpu-closed-dp-component',
           'native-cpu-closed-dp-provider', 'native-cpu-closed-dp-cohort']
FAILURES = {'out/native-cpu-canonical-closed-dp/value/owner.json',
            'out/native-cpu-closed-dp-component/build/owner.json',
            'out/native-cpu-closed-dp-cohort/control-8/owner.json'}


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not A.exists()
owners, sources, products = {}, [Path(__file__).resolve()], []
for folder in FOLDERS:
    for path in (R / 'out' / folder).rglob('*'):
        if not path.is_file() or '__pycache__' in path.parts:
            continue
        if path.suffix in ('.so', '.o', '.bc'):
            products.append(path)
        else:
            sources.append(path)
        if path.name == 'owner.json':
            name, row = str(path.relative_to(R)), doc(path)
            assert row['status'] == 'terminal' and row['passed'] is (name not in FAILURES)
            assert row['child_exit'] == (1 if name in FAILURES else 0)
            assert row['source_unchanged'] and row['whole_inputs_unchanged']
            assert not row['remaining_owned_processes'] and not row['cleanup_errors']
            assert row['observed_births'] and all(item['absent'] for item in row['observed_births'])
            assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()
            owners[name] = dict(identity=pin(path), passed=row['passed'], source=row['source'],
                                wall_seconds=row['wall_seconds'], peak_sampled_rss_kib=row['peak_sampled_rss_kib'])
assert len(owners) == 16 and {name for name, row in owners.items() if not row['passed']} == FAILURES
observation_path = R / 'out/native-cpu-closed-dp-cohort/control-8/native-observation.json'
observation = doc(observation_path)
assert observation['failures'] == ['shared-value entry: AssertionError()']
assert observation['inferior_exit'] is None
assert not (R / 'out/native-cpu-closed-dp-cohort/control-8-acceptance.json').exists()
last = observation['optional_math_calls'][-1]
assert last['function'] == 'sirius_local_value_dense_add8' and not last['returned']
assert last['raw_caller']['name'] == 'sirius_closed_dp_packet'
proof_path = R / 'out/native-cpu-closed-dp-provider/provider/linked-native-proof.json'
proof = doc(proof_path)
assert proof['passed'] and proof['closed_dp_interface']['passed']
library = R / proof['library']
assert pin(library) == proof['library_identity']
products.append(library)
for name in ['canonical_private_dp.py', 'canonical_private_dp_contract.json',
             'derive_closed_dp_passes.py', 'constrain_closed_dp.cpp', 'closed_dp_scope_contract.json',
             'closed_dp_value_abi.cpp', 'closed_dp_value_contract.json',
             'derive_closed_dp_context.py', 'closed_dp_context.cpp', 'closed_dp_context_contract.json']:
    sources.append(V / name)
base = R / 'bin/linux-gcc/native-cpu-canonical-closed-dp'
sources += [path for path in base.rglob('*') if path.is_file() and path.suffix in ('.ll', '.ispc')]
products += [path for path in base.iterdir() if path.is_file() and path.suffix in ('.bc', '.o', '.so')]
products += [R / 'bin/linux-gcc/native-cpu-local-effects' / name for name in ('math.bc', 'math.o')]
products += [R / 'bin/linux-gcc/native-cpu-closed-dp-component' / name for name in ('adapter-only.o', 'rows_adapter.o')]
predecessors = {}
for name in ['attestations/qualification/844dbeda/native-cpu-shared-value-original-frame-attempt/MANIFEST.json',
             'attestations/qualification/359307b1/native-cpu-register-value-compiler-gate/MANIFEST.json',
             'attestations/qualification/b2cbd4db/native-cpu-private-context-interface-gates/MANIFEST.json']:
    predecessors[name] = pin(R / name)
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
manifest = dict(status='CLOSED_DP_COMPILER_PASS_C8_OBSERVER_STOP_NO_NUMERICAL_VERDICT',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    candidate_revision='4332fb950b3ced25bd28d5792bb4c478b0dc78b3', owners=owners,
    retained_receipt_failures=sorted(FAILURES), predecessors=predecessors, bundles=bundles,
    compiler_context_gate=doc(R / 'out/native-cpu-canonical-closed-dp/cleanup-result.json'),
    original_C8_observation=pin(observation_path),
    original_C8_failures=observation['failures'],
    numerical_verdict=False,
    scope='One original C8 attempt for the new closed canonical SSA DP/local9 strict cc92 math stopped at its first local dense_add entry because the passive observer could not unwind beyond the correct closed DP caller to prove the outer parent/job frames. The inferior did not complete and no32-digest numerical verdict exists; this is a failed required observation gate. AU T/E/shared9 and all public/fallback branches remain exact. The120832-byte private Context is gone from optimized DP IR; one small alloca/five ordinary loads/zero ordinary stores remain. Original constrained arithmetic, masks, tables, operation/root order and FP quarantine are preserved. Both successful compiler/link products whose receipt readers failed were accepted through distinct read-only audits with original failures preserved, without repeating compilation/linking. No unchanged workload retry, ABBA22, Decimal, all8, original frame, physical hardware, scientific or release qualification. Ultimate goal and all eight subgoals remain open.')
with (A / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)), owners=len(owners),
    passed=sum(row['passed'] for row in owners.values()), preserved_failed_owners=len(FAILURES),
    bundles={key:value['bytes'] for key,value in bundles.items()})))
