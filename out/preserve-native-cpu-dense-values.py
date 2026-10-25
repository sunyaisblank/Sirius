"""Preserve the terminal R candidate once, with complete member hashes and modes."""
from pathlib import Path
import hashlib
import json
import os
import re
import stat
import subprocess
import zipfile

R = Path.cwd()
W = R / 'out/native-cpu-dense-values'
C = R / 'out/native-cpu-dense-values-cohort'
D = R / 'out/native-cpu-dense-values-compatibility'
A = R / 'out/native-cpu-dense-arithmetic'
Q = R / 'out/native-cpu-direct-scratch-cohort'


def doc(path):
    return json.loads(path.read_text())


def identity(path):
    before = path.stat()
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns, before.st_mode) == (after.st_size, after.st_mtime_ns, after.st_mode)
    return {'bytes': after.st_size, 'sha256': digest,
            'mode': stat.S_IMODE(after.st_mode), 'mtime_ns': after.st_mtime_ns}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert doc(C / 'control-8-acceptance.json')['passed']
assert doc(Q / 'benefit/readers/benefit.json')['all22_numeric_pass']
assert doc(Q / 'benefit/readers/decimal-agreement.json')['all22_pass']
assert doc(Q / 'controls/readers/controls-1-8-acceptance.json')['all8_original_controls_pass']
assert doc(D / 'result.json')['passed']
frame = doc(C / 'cold-frame/result.json')
frame_owner = doc(C / 'cold-frame/frame/owner.json')
assert frame_owner['status'] == 'terminal'
scope = 'accepted-original-frame' if frame_owner['passed'] else 'frame-timeout'
assert frame_owner['passed'] or frame_owner['stop_reason'] == 'diagnostic_timeout'
destination = R / 'attestations/qualification' / revision[:8] / ('native-cpu-dense-values-' + scope)
assert not destination.exists()
owners = {}
boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
for root in (W, C, D):
    for path in sorted(root.rglob('owner.json')):
        owner = doc(path)
        assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        if not owner['passed']:
            assert path == C / 'cold-frame/frame/owner.json'
        for row in owner['observed_births']:
            process = Path('/proc') / str(row['pid']) / 'stat'
            if boot == owner['boot_id'] and process.exists():
                try:
                    actual = process.read_text().rsplit(')', 1)[1].split()[19]
                except FileNotFoundError:
                    continue
                assert actual != str(row['saved_start_ticks']), (path, row)
        owners[str(path.relative_to(R))] = {key: owner.get(key) for key in
            ('source', 'status', 'passed', 'child_exit', 'wall_seconds', 'peak_sampled_rss_kib', 'stop_reason')}
assert len(owners) == 8
assert frame_owner['source']['revision'] == revision
progress = re.findall(rb'Tile\s+(\d+)\s*/\s*(\d+)', (C / 'cold-frame/frame/stdout.log').read_bytes())
compiled = doc(W / 'ispc-values-compile-result.json')
provider = doc(W / 'provider/provider-build-result.json')
diagnostic = doc(D / 'link-result.json')
products = [R / row['object']['path'] for row in compiled['objects']]
products += [R / compiled['merged_adapter']['path']]
bridge = doc(W / 'bridge-build-result.json')
products += [R / bridge['object']['path']]
products += [R / provider[key] for key in ('library', 'ICD')]
products += [R / diagnostic['executable']['path']]
assert len(products) == len(set(products)) == 8
source_files = {str(path.relative_to(R)): path for root in (W, C, D)
                for path in root.rglob('*') if path.is_file()}
for path in (R / 'bin/linux-gcc/native-cpu-dense-values/ispc-source').glob('*.ispc'):
    source_files[str(path.relative_to(R))] = path
for path in (Path(__file__).resolve(),):
    source_files[str(path.relative_to(R))] = path
assert len([name for name in source_files if name.endswith('.ispc')]) == 3
product_files = {str(path.relative_to(R)): path for path in products}
runtime = doc(C / 'runtime-transition.json')
assert {key: identity(Path(runtime['path']))[key] for key in ('bytes', 'sha256')} == runtime['current_runtime']
destination.mkdir(parents=True)


def archive(name, files):
    path = destination / name
    members = {}
    with zipfile.ZipFile(path, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as output:
        for key, source in sorted(files.items()):
            before = identity(source)
            info = zipfile.ZipInfo.from_file(source, 'files/' + key)
            info.compress_type = zipfile.ZIP_DEFLATED
            output.writestr(info, source.read_bytes())
            assert identity(source) == before
            members['files/' + key] = {'source': str(source), **before}
    with zipfile.ZipFile(path) as output:
        assert output.testzip() is None and set(output.namelist()) == set(members)
        for key, expected in members.items():
            data = output.read(key)
            assert len(data) == expected['bytes'] and hashlib.sha256(data).hexdigest() == expected['sha256']
            assert stat.S_IMODE(output.getinfo(key).external_attr >> 16) == expected['mode']
            assert identity(Path(expected['source'])) == {k: expected[k] for k in ('bytes', 'sha256', 'mode', 'mtime_ns')}
    return {'file': name, 'identity': identity(path), 'members': members}


bundles = [archive('source-build-and-trial.zip', source_files), archive('actual-products.zip', product_files)]
references = {}
for name in (
    'attestations/qualification/2d971dff/native-cpu-nir-dead-progress-frame-timeout/MANIFEST.json',
    'attestations/qualification/37cd12be/native-cpu-upload-aware-staging-rejected-transport12/MANIFEST.json',
    'attestations/qualification/1fdaa56a/native-cpu-worker-occupancy-frame-timeout/MANIFEST.json',
    'attestations/qualification/603a155e/native-cpu-accumulate-inline-frame-timeout/MANIFEST.json',
    'attestations/qualification/4595f5a6/native-cpu-product-inline-frame-timeout/MANIFEST.json',
    'attestations/qualification/bee63496/native-cpu-task-cardinality-frame-timeout/MANIFEST.json',
    'attestations/qualification/7a8e6d48/native-cpu-direct-scratch-frame-timeout/MANIFEST.json',
    'attestations/qualification/6c79517c/native-cpu-static-registers-frame-timeout/MANIFEST.json',
    'attestations/qualification/fddb9175/native-cpu-predicated-products-frame-timeout/MANIFEST.json',
    'attestations/qualification/2e8224da/native-cpu-dense-arithmetic-frame-timeout/MANIFEST.json',
):
    references[name] = identity(R / name)
manifest = {'schema': 'sirius-native-cpu-dense-values-terminal-v1',
            'qualifying_source_revision': revision, 'accepted': bool(frame['accepted']), 'owners': owners, 'archives': bundles,
            'additional_qualification_claimed': False,
            'scope': 'AD routes remaining compound Add/Multiply through initialized explicit word buffers and packing-only C++ adapters into exact checked AC arithmetic. Original AC math/adapter object and all its function opcodes are unchanged; no arithmetic compilation/source change. AC Store calls remain, with all9380tuples and other original row/control/refusal/rollback/publication. Masks and operand/output snapshots preserve aliases and inactive words. New baseline packing-only object and actual raw-to-value-to-original-math target graph joined. Fresh full38/all187960words and original C8/all32+51markers+54joins passed; actual original frame result is decisive. Broader Z ABBA22+Decimal/all8 remain predecessor evidence. No independent math/science/release/main/780M/large-scene/original-July completion claim.',
            'original_frame': frame, 'runtime_transition': doc(C / 'runtime-transition.json'),
            'last_reported_tiles': list(map(int, progress[-1])) if progress else None,
            'additional_passive_profile_executed': False,
            'profile_scope': 'No profiling study for AD. Historical R profile supports the arithmetic premise only; it is not current Z execution or performance evidence.',
            'dense_CPP_ABI_scope': {'focused_owner': identity(A / 'cpp/run/owner.json'), 'focused_output': doc(A / 'cpp/run/stdout.log'), 'CXX_build': identity(A / 'cpp-build-result.json'), 'CXX_code': identity(A / 'cpp-code-inspection.json')},
            'diagnostic_build_scope': 'Zero new CPP compilations or baseline objcopy operations. Reused exact verified R driver and renamed Q baseline; one new link and one fresh full38/all187960word comparison passed. Historical R accounting failures remain historical.',
            'predecessor_Z_broader_evidence': {name: identity(Q / name) for name in ('benefit/readers/benefit.json', 'benefit/readers/decimal-agreement.json', 'controls/readers/controls-1-8-acceptance.json')},
            'legacy_human_scope_labels': 'Inherited comments may retain old stage labels; frozen paths/object hashes/native PCs/module joins identify actual Z. Original ABBA/frame scopes identify Z. No new AD broad-suite execution; actual current paths/hashes and fresh C8 identify AD. Z broader results are preserved separately.',
            'new_product_members': 8, 'current_system_zlib_recovery_member': False,
            'reused_driver_baseline_runtime_recovery': 'R603a155e actual-products.zip preserves the exact driver object, renamed baseline and current f61677c7 system zlib; no duplicate copies here.',
            'reused_native_objects': doc(W / 'provider/source-contract.json')['native_objects'][:4],
            'recovery_manifests': references,
            'resource_payloads': 'Original32 cache assets and their existing recovery bundle reused; frame aliases can be retired only after this archive is independently verified and committed/pushed.'}
with (destination / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps({'bundle': str(destination.relative_to(R)), 'owners': len(owners),
                  'archives': [{k: v for k, v in item.items() if k != 'members'} for item in bundles]}))
