"""Preserve the terminal AJ candidate once, with complete member hashes and modes."""
from pathlib import Path
import hashlib
import json
import os
import re
import stat
import subprocess
import zipfile

R = Path.cwd()
W = R / 'out/native-cpu-normalization-tags'
C = R / 'out/native-cpu-normalization-tags-cohort'
D = R / 'out/native-cpu-normalization-tags-compatibility'
A = W
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
assert doc(C / 'affected-control/result.json')['passed']
assert doc(Q / 'benefit/readers/benefit.json')['all22_numeric_pass']
assert doc(Q / 'benefit/readers/decimal-agreement.json')['all22_pass']
assert doc(Q / 'controls/readers/controls-1-8-acceptance.json')['all8_original_controls_pass']
assert doc(D / 'result.json')['passed']
assert doc(W/'predicate/result.json')['passed']
frame = doc(C / 'cold-frame/result.json')
frame_owner = doc(C / 'cold-frame/frame/owner.json')
assert frame_owner['status'] == 'terminal'
scope = 'accepted-original-frame' if frame_owner['passed'] else 'frame-timeout'
assert frame_owner['passed'] or frame_owner['stop_reason'] == 'diagnostic_timeout'
destination = R / 'attestations/qualification' / revision[:8] / ('native-cpu-normalization-tags-' + scope)
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
            if path == W/'predicate/link/owner.json':
                admitted=doc(W/'predicate/link-result.json')
                assert admitted['passed'] and admitted['original_link_owner_passed'] is False
                assert {k:identity(path)[k] for k in ('bytes','sha256')} == admitted['original_link_owner']
                assert owner['child_exit']==1 and owner['stop_reason'] is None
            else:
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
assert len(owners) == 12
assert frame_owner['source']['revision'] == revision
progress = re.findall(rb'Tile\s+(\d+)\s*/\s*(\d+)', (C / 'cold-frame/frame/stdout.log').read_bytes())
compiled = doc(W / 'ispc-dense-compile-result.json')
provider = doc(W / 'provider/provider-build-result.json')
diagnostic = doc(D / 'link-result.json')
products = [R / row['object']['path'] for row in compiled['objects']]
products += [R / compiled['merged_adapter']['path']]
cpp = doc(W / 'cpp-build-result.json')
products += [R / row['object']['path'] for row in cpp['compilations']]
products += [R / cpp['executable']['path']]
products += [R / provider[key] for key in ('library', 'ICD')]
products += [R / diagnostic['executable']['path']]
products += [R/compiled['predicate_compilation']['object']['path'],R/doc(W/'predicate/link-result.json')['executable']['path']]
assert len(products) == len(set(products)) == 12
source_files = {str(path.relative_to(R)): path for root in (W, C, D)
                for path in root.rglob('*') if path.is_file()}
for path in (R / 'bin/linux-gcc/native-cpu-normalization-tags/ispc-source').glob('*.ispc'):
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
    'attestations/qualification/4e49d817/native-cpu-dense-values-frame-timeout/MANIFEST.json',
    'attestations/qualification/46bfaaf8/native-cpu-original-frame-route/MANIFEST.json',
    'attestations/qualification/14b5567d/native-cpu-masked-arithmetic-frame-timeout/MANIFEST.json',
    'attestations/qualification/6176354b/native-cpu-demand-submit-frame-timeout/MANIFEST.json',
    'attestations/qualification/8657e016/native-cpu-arithmetic-closure-frame-timeout/MANIFEST.json',
    'attestations/qualification/9275d0aa/native-cpu-direct-compound-frame-timeout/MANIFEST.json',
    'attestations/qualification/f878a18f/native-cpu-normalized-arithmetic-frame-timeout/MANIFEST.json',
):
    references[name] = identity(R / name)
manifest = {'schema': 'sirius-native-cpu-normalization-tags-terminal-v1',
            'qualifying_source_revision': revision, 'accepted': bool(frame['accepted']), 'owners': owners, 'archives': bundles,
            'additional_qualification_claimed': False,
            'scope': 'AJ caches the unchanged normalization predicate for the same complete private register value. Whole context initialization covers tags; both generic/uniform partial writes invalidate them and complete writes/direct stores refresh them after payload writes. No assumed normalization theorem, extra published limb, persistent/input/shader-cache change, stronger eligibility, changed arithmetic/formula/FP/ISA/precision/requests/fields/RHS/trials/quartics/order/refusal/rollback/governors. All9380ordered literal calls and public layouts remain. Original AI CPP source is an exact prefix. Fresh370440call44452800payload-word8890560tag-word original-scalar/independent-nextafter fixture passes, including aliases/masks/false-negative tags. The exact shared ISPC predicate passes1764calls42336tagwords/211680unchanged-inputwords. Initial object reader incorrectly required tagged Sqrt in DP with no static Sqrt sites; original frozen reader retained and corrected reader applies actual canonical sites, no compiler repeat. Predicate link completed then a reader selected driver.o instead of math.o; corrected association reuses the same ELF and original failed owner remains failed. No link/compiler/workload repeated for that correction. Exact H host/frame reused; current full38/C8/affected control5 and sole original frame decide their finite scopes. Z all8/ABBA22+Decimal remain predecessor evidence; no full scientific/Mandatory/release acceptance.',
            'original_frame': frame, 'runtime_transition': doc(C / 'runtime-transition.json'),
            'last_reported_tiles': list(map(int, progress[-1])) if progress else None,
            'additional_passive_profile_executed': doc(C/'cold-frame/profile-controller.json')['observer_started'],
            'current_profile_controller':doc(C/'cold-frame/profile-controller.json'),
            'profile_scope': 'One predeclared fixed3s existing USER-PC capture during the same original frame; no extra render. Original frame wall includes observer interference, with observer footprint recorded separately and not subtracted. Samples are a bounded window, not uninstrumented speed or whole-frame/kernel/off-CPU causality.',
            'dense_CPP_ABI_scope': {'focused_owner': identity(A / 'cpp/run/owner.json'), 'focused_output': doc(A / 'cpp/run/stdout.log'), 'CXX_build': identity(A / 'cpp-build-result.json'), 'CXX_code': identity(A / 'cpp-code-inspection.json')},
            'reused_host_build': identity(R/'out/native-cpu-owned-lower-increments/host-build-result.json'),
            'fresh_affected_control5': doc(C/'affected-control/result.json'),
            'diagnostic_build_scope': 'Zero new CPP compilations or baseline objcopy operations. Reused exact verified R driver and renamed Q baseline; one new link and one fresh full38/all187960word comparison passed. Historical R accounting failures remain historical.',
            'predecessor_Z_broader_evidence': {name: identity(Q / name) for name in ('benefit/readers/benefit.json', 'benefit/readers/decimal-agreement.json', 'controls/readers/controls-1-8-acceptance.json')},
            'legacy_human_scope_labels': 'Some inherited comments retain predecessor stage names; exact frozen paths/hashes/native PCs/module joins identify AJ. Fresh focused/predicate/ISPC/full38/C8/affected-control5/frame owners are recorded here. Z broader results remain separately preserved predecessor evidence.',
            'shared_predicate_result':doc(W/'predicate/result.json'),'predicate_link_admission':doc(W/'predicate/link-result.json'),
            'new_product_members': 12, 'current_system_zlib_recovery_member': False,
            'reused_driver_baseline_runtime_recovery': 'R603a155e actual-products.zip preserves the exact driver object, renamed baseline and current f61677c7 system zlib; no duplicate copies here.',
            'reused_native_objects': doc(W / 'provider/source-contract.json')['native_objects'][:4],
            'reused_H_host_recovery': 'AH f878a18f actual-products.zip preserves the exact H host archive/executables/integrator object; its source ZIP preserves the H source/build/ownership proof. No duplicate H products here.',
            'recovery_manifests': references,
            'resource_payloads': 'Original32 cache assets and their existing recovery bundle reused; frame aliases can be retired only after this archive is independently verified and committed/pushed.'}
with (destination / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps({'bundle': str(destination.relative_to(R)), 'owners': len(owners),
                  'archives': [{k: v for k, v in item.items() if k != 'members'} for item in bundles]}))
