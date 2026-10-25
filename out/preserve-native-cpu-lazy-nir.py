"""Preserve V's single original frame rejection and its preceding gates."""
from pathlib import Path
import ast
import hashlib
import json
import re
import stat
import subprocess
import zipfile

R = Path.cwd()
W = R / 'out/native-cpu-lazy-nir'
C = R / 'out/native-cpu-lazy-nir-cohort'
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
basis = R / 'out/preserve-native-cpu-zero-term.py'
tree = ast.parse(basis.read_text())
functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in ('doc', 'identity', 'archive')]
assert len(functions) == 3
exec(compile(ast.Module(body=functions, type_ignores=[]), str(basis), 'exec'))
assert doc(C / 'control-8-acceptance.json')['passed']
assert doc(C / 'benefit/readers/benefit.json')['all22_numeric_pass']
assert doc(C / 'benefit/readers/decimal-agreement.json')['all22_pass']
assert doc(C / 'controls/readers/controls-1-8-acceptance.json')['all8_original_controls_pass']
frame = doc(C / 'cold-frame/result.json')
failed = C / 'cold-frame/frame/owner.json'
timeout = doc(failed)
assert timeout['status'] == 'terminal' and not timeout['passed']
assert timeout['stop_reason'] == 'diagnostic_timeout' and timeout['source']['revision'] == revision
owners = {}
boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
for root in (W, C):
    for path in sorted(root.rglob('owner.json')):
        owner = doc(path)
        assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['passed'] or path == failed
        for row in owner['observed_births']:
            process = Path('/proc') / str(row['pid']) / 'stat'
            if boot == owner['boot_id'] and process.exists():
                try:
                    actual = process.read_text().rsplit(')', 1)[1].split()[19]
                except FileNotFoundError:
                    continue
                assert actual != str(row['saved_start_ticks'])
        owners[str(path.relative_to(R))] = {key: owner.get(key) for key in
            ('source', 'status', 'passed', 'child_exit', 'wall_seconds', 'peak_sampled_rss_kib', 'stop_reason')}
assert len(owners) == 15
progress = re.findall(rb'Tile\s+(\d+)\s*/\s*(\d+)', (C / 'cold-frame/frame/stdout.log').read_bytes())
destination = R / 'attestations/qualification' / revision[:8] / 'native-cpu-lazy-nir-frame-timeout'
assert not destination.exists()
destination.mkdir(parents=True)
sources = {str(path.relative_to(R)): path for root in (W, C) for path in root.rglob('*') if path.is_file()}
for path in (Path(__file__).resolve(), basis):
    sources[str(path.relative_to(R))] = path
mesa = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
for name in ('src/gallium/drivers/llvmpipe/lp_screen.c', 'src/gallium/drivers/llvmpipe/lp_state_cs.c',
             'src/gallium/drivers/llvmpipe/lp_state_cs.h', 'src/gallium/frontends/lavapipe/lvp_device.c'):
    path = mesa / name
    sources[str(path.relative_to(R))] = path
provider = doc(W / 'provider/provider-build-result.json')
products = {str((R / provider[key]).relative_to(R)): R / provider[key] for key in ('library', 'ICD')}
assert len(products) == 2
bundles = [archive('source-build-and-trial.zip', sources), archive('actual-products.zip', products)]
manifest = {
    'schema': 'sirius-native-cpu-lazy-nir-terminal-v1', 'qualifying_source_revision': revision,
    'accepted': False, 'owners': owners, 'archives': bundles, 'original_frame': frame,
    'last_reported_tiles': list(map(int, progress[-1])), 'new_product_members': 2,
    'scope': 'V defers native marked compute NIR optimization to first real variant generation, before IR cache serialization; other shaders retain eager finalization. Actual U scalar/S row arithmetic, Q scheduler, host, fixture and original scene unchanged. Incremental provider build and full898function/module joins passed, reused exact U full38 checker joined, original C8/all32/51markers/54joins, sole ABBA22 plus independent Decimal, and all8 original controls passed. The sole original128x128 one-SPP software Vulkan frame timed out under unchanged300s/4GiB/2GiB governors. No unchanged retry, full science/release/main or adoption claim.',
    'control3_observer_limit': 'The first64 finalizer entries were all graphics shader_stage4/native_stage0 and exhausted the capture cap. The observer recorded the first EndpointRows return true, but does not demonstrate compute tag classification, lazy skip, fallback, or the full timeout cause. Original C3 assertions passed; no repeated capture.',
    'additional_frame_profile_executed': False,
    'native_objects': doc(W / 'provider/source-contract.json')['native_objects'],
    'reused_U_products_recovery': 'attestations/qualification/fb50bb2a/native-cpu-scalar-zero-control3-timeout/MANIFEST.json',
    'reused_S_row_products_recovery': 'attestations/qualification/115b86f3/native-cpu-zero-term-frame-timeout/MANIFEST.json',
    'reused_driver_baseline_runtime_recovery': 'attestations/qualification/603a155e/native-cpu-accumulate-inline-frame-timeout/MANIFEST.json',
    'runtime_transition': doc(C / 'runtime-transition.json'),
    'resource_payloads': 'Original32 cache assets/recovery preserved. Retire only restored hardlink aliases after independent bundle verification and commit/push.'}
with (destination / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps({'bundle': str(destination.relative_to(R)), 'owners': len(owners),
    'archives': [{'file': row['file'], 'bytes': row['identity']['bytes'], 'sha256': row['identity']['sha256'],
                  'members': len(row['members'])} for row in bundles]}))
