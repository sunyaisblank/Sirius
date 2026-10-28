"""Seal the one coupled shared-value attempt after its first terminal gate."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R = Path.cwd()
V = R / 'out/native-cpu-strict-ir-linking'
C = R / 'out/native-cpu-shared-value-component'
P = R / 'out/native-cpu-shared-value-provider'
Q = R / 'out/native-cpu-shared-value-cohort'
A = R / 'attestations/qualification/844dbeda/native-cpu-shared-value-original-frame-attempt'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not A.exists()
structural = [V / name for name in ('shared-plugin', 'shared-math-lower', 'shared-math-opt',
    'shared-math-native', 'shared-dp-lower', 'shared-dp-opt', 'shared-dp-native',
    'shared-transport-lower', 'shared-transport-opt', 'shared-transport-native',
    'shared-endpoint-lower', 'shared-endpoint-opt', 'shared-endpoint-cleanup-opt',
    'shared-endpoint-cleanup-native')]
owners = structural + [C / 'build', P / 'provider/configure', P / 'provider/build', Q / 'control-8']
owners += [Q / 'benefit/runtime' / name for name in ('a1', 'b1', 'b2', 'a2')]
owners += [Q / 'controls/runtime' / ('control-' + str(i)) for i in range(1, 8)]
owners += [Q / 'cold-frame/frame']
history = {}
for folder in owners:
    path = folder / 'owner.json'
    owner = json.loads(path.read_text())
    assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    for birth in owner['observed_births']:
        try:
            current = (Path('/proc') / str(birth['pid']) / 'stat').read_text().rsplit(')', 1)[1].split()[19]
        except FileNotFoundError:
            continue
        assert current != birth['saved_start_ticks'], (folder, birth)
    history[str(folder.relative_to(R))] = dict(path=str(path.relative_to(R)), **pin(path),
        passed=owner['passed'], source=owner['source'], seconds=owner['wall_seconds'],
        peak_rss_kib=owner['peak_sampled_rss_kib'], stop_reason=owner['stop_reason'])
    if not owner['passed']:
        assert folder in (V / 'shared-endpoint-opt', Q / 'cold-frame/frame')
assert len(history) == 30
source_paths = [Path(__file__).resolve()]
source_paths += [p for parent in (C, P, Q) for p in parent.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
source_paths += [p for p in V.iterdir() if p.is_file() and
                 (p.name.startswith(('shared', 'endpoint_cleanup')) or p.name in
                  ('scope_contract.json', 'remaining_scope_review.json', 'toolchain.json'))]
source_paths += [p for folder in structural for p in folder.iterdir() if p.is_file()]
source_paths += [R / 'out/native-cpu-exact-product-zero-residual/retained8.cpp',
                 R / 'out/native-cpu-inline-private-bank/packinghelper.isph']
source_paths += [R / ('bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/' + s + '.ispc')
                 for s in ('transport', 'endpoint', 'dopri_phase')]
products = [R / json.loads((P / 'provider/provider-build-result.json').read_text())['library']]
products += [R / json.loads((C / 'build-result.json').read_text())['output']['path']]
products += [R / 'bin/linux-gcc/native-cpu-shared-value-component/adapter-only.o']
for name in ('shared-dp', 'shared-transport', 'shared-math', 'shared-endpoint-cleanup'):
    products += [R / json.loads((V / (name + '-native-result.json')).read_text())['output']['path']]
    products += [R / json.loads((V / (name + '-opt-result.json')).read_text())['output']['path']]
source_paths = [path for path in source_paths if path not in set(products)]
A.mkdir(parents=True)
bundles = {}
for label, paths in (('source-and-evidence', source_paths), ('products', products)):
    entries = {str(p.relative_to(R)): pin(p) for p in sorted(set(paths))}
    target = A / (label + '.zip')
    with zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zipped:
        for name in entries:
            zipped.write(R / name, name)
    with zipfile.ZipFile(target) as zipped:
        assert set(zipped.namelist()) == set(entries)
        for name, expected in entries.items():
            with zipped.open(name) as stream:
                assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected['sha256'], name
            assert zipped.getinfo(name).file_size == expected['bytes']
            assert pin(R / name) == expected, name
    assert target.stat().st_size < 100 * 1024 * 1024
    bundles[label] = dict(path=str(target.relative_to(R)), **pin(target), members=entries)
frame = json.loads((Q / 'cold-frame/frame/owner.json').read_text())
manifest = dict(status='ORIGINAL_FRAME_PASS' if frame['passed'] else 'ORIGINAL_FRAME_FAILED_STOP',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    execution_revision='844dbeda14ab6fb3ef09fddf0c98f5b254bb89c2', owners=history,
    phase_passes=dict(Counter('passed' if x['passed'] else 'failed' for x in history.values())),
    bundles=bundles, all3_private_payload_banks_eliminated=8241, all3_private_tag_banks_eliminated=8181,
    shared_core_count=9, static_native_call_joins=8241,
    source_authority='Exact owner.source commits and per-action frozen inputs remain authoritative; no old proof/receipt is rewritten. Endpoint full second O2 timeout is preserved; focused cleanup of already-ISPC-O2 input passed without repeating DP/T or raising limits.',
    component_interface=json.loads((C / 'abi-review.json').read_text()),
    C8=json.loads((Q / 'control-8-acceptance.json').read_text()),
    benefit=json.loads((Q / 'benefit/readers/benefit.json').read_text()),
    decimal=json.loads((Q / 'benefit/readers/decimal-agreement.json').read_text()),
    all8=json.loads((Q / 'controls/readers/controls-1-8-acceptance.json').read_text()),
    final_frame_result=json.loads((Q / 'cold-frame/result.json').read_text()),
    final_frame_owner=frame,
    disposable_cleanup=json.loads((Q / 'cold-frame/disposable-cleanup.json').read_text()),
    reader_correction='The original inherited reader generated an AT evidence root. Separate read_terminal.py corrects only that literal to AU; all attempted workload inputs, owner receipts and original reader assertions remain unchanged. No workload rerun.',
    prior_strict_conversion_gate=dict(path='attestations/qualification/9cdafa03/native-cpu-strict-value-interface-compiler-gate/MANIFEST.json',
        **pin(R / 'attestations/qualification/9cdafa03/native-cpu-strict-value-interface-compiler-gate/MANIFEST.json')),
    scope='One materially changed shared-value implementation, same original C8/all32, soleABBA22 plus independentDecimal, seven original controls/all8, and single original128x128/oneSPP frame. No changed scenes/counts/layouts/allocations/precision/thresholds/work limits/cache/governors, hidden warmups/retries or scientific/hardware-domain/release acceptance transfer.')
write(A / 'MANIFEST.json', manifest)
print(json.dumps(dict(sealed=str(A.relative_to(R)), status=manifest['status'],
    owners=len(history), phase_passes=manifest['phase_passes'], bundles={k: v['bytes'] for k, v in bundles.items()})))
