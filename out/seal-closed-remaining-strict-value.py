"""Preserve the first strict/value compiler artifacts; no numerical transfer."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R = Path.cwd()
W = R / 'out/native-cpu-canonical-closed-remaining'
V = R / 'out/native-cpu-strict-ir-linking'
A = R / 'attestations/qualification/0de5e321/native-cpu-closed-remaining-strict-value-gates'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not A.exists()
sources = [Path(__file__).resolve()]
sources += [path for path in W.iterdir() if path.is_file() and path.suffix in ('.py', '.json')]
sources += [V / 'derive_closed_remaining_passes.py', V / 'private_finite_closure.json']
products, owners, results = [], {}, {}
for stage in ('transport', 'endpoint'):
    B = R / ('bin/linux-gcc/native-cpu-canonical-closed-' + stage)
    sources += [V / ('constrain_closed_' + stage + '.cpp'), V / ('closed_' + stage + '_scope_contract.json'),
        V / ('closed_' + stage + '_value_abi.cpp'), V / ('closed_' + stage + '_value_contract.json')]
    sources += [B / (stage + '.' + mode + '.ll') for mode in ('strict', 'value')]
    products += [B / name for name in ('scope.so', 'value.so', stage + '.strict.bc', stage + '.value.bc')]
    for mode in ('scope', 'value'):
        for variant in ('-plugin', ''):
            action = mode + variant + '-' + stage
            folder = W / action
            row = doc(folder / 'owner.json')
            assert row['status'] == 'terminal' and row['passed'] and row['child_exit'] == 0
            assert row['source_unchanged'] and row['whole_inputs_unchanged']
            assert not row['remaining_owned_processes'] and not row['cleanup_errors']
            assert row['observed_births'] and all(item['absent'] for item in row['observed_births'])
            assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
            sources += [path for path in folder.iterdir() if path.is_file()]
            owners[action] = dict(path=str((folder / 'owner.json').relative_to(R)), identity=pin(folder / 'owner.json'),
                source=row['source'], wall_seconds=row['wall_seconds'], peak_sampled_rss_kib=row['peak_sampled_rss_kib'])
            receipt = W / (stage + '-' + mode + variant + '-result.json')
            result = doc(receipt)
            assert pin(R / result['output']['path']) == {key:result['output'][key] for key in ('bytes', 'sha256')}
            results[action] = result
assert len(owners) == 8
predecessors = {}
for name in ['attestations/qualification/d970f13e/native-cpu-closed-remaining-first-emit/MANIFEST.json',
    'attestations/qualification/620a0c6e/native-cpu-closed-dp-unwind-original-c8/MANIFEST.json',
    'attestations/qualification/620a0c6e/native-cpu-closed-dp-unwind-original-c8/ALLOCATION_CORRECTION.json']:
    predecessors[name] = pin(R / name)
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
manifest = dict(status='STRICT_VALUE_COMPILER_GATES_PASS_CONTEXT_AND_NUMERICAL_GATES_PENDING',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    predecessors=predecessors, owners=owners, results=results, bundles=bundles,
    scope='First builds/executions of exact Transport/Endpoint strict and value passes, eight PASS compiler owners. Both full modules/source/layouts/casts/caller/table/FP inventories are pinned. T6 full8 casts are guarded as defined exactly representable integers; E original scale/comparison checks remain. Original quiet/signaling, dynamic rounding, strict exceptions, fneg/fabs and all2698/3954 raw call sites remain accounted for; pre-inline typed wrapper calls are8/9. Same local9 math definitions are reused without recompilation, with cc92/inaccessibleRW and unconditional private ordinary-memory closure. Original transaction staging remains mandatory: T150368B and E122528B; DP132768B staging is separately corrected in the linked append-only record. No context promotion/native/FDE/provider/runtime/numerical/benefit/all8/ABBA/frame/scientific/hardware/release acceptance is inferred. All eight ultimate-goal subgoals remain open.')
with (A / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)), passed_owners=8,
    bundles={name:row['bytes'] for name,row in bundles.items()})))
