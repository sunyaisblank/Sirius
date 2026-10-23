"""Build and stage the private native CPU bridge once; do not load the provider."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
P = W / 'provider'
B = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
destination = R / 'bin/linux-gcc/native-cpu-slang-bridge/provider'


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


plan = json.loads((P / 'incremental-plan.json').read_text())
assert not destination.exists()
for name, identity in plan['reused_build_inputs'].items():
    assert pin(R / name) == identity, name
target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'
assert pin(target) == plan['previous_target']
subprocess.run(['ninja', '-C', str(B), '-j4',
                'src/gallium/targets/lavapipe/libvulkan_lvp.so'], check=True)
for name, identity in plan['reused_build_inputs'].items():
    assert pin(R / name) == identity, name
destination.mkdir(parents=True)
library = destination / 'libvulkan_lvp.so'
shutil.copy2(target, library)
assert pin(library) == pin(target) and pin(target) != plan['previous_target']
icd = destination / 'lvp_icd.json'
icd.write_text(json.dumps({'file_format_version': '1.0.0',
    'ICD': {'library_path': str(library), 'api_version': '1.4.305'}}, indent=2) + '\n')
(P / 'provider-build-result.json').write_text(json.dumps({
    'passed': True, 'library': str(library.relative_to(R)), 'library_identity': pin(library),
    'ICD': str(icd.relative_to(R)), 'ICD_identity': pin(icd),
    'actual_build_target': str(target.relative_to(R)), 'target_identity': pin(target),
    'reused_build_inputs_unchanged': len(plan['reused_build_inputs']),
    'native_objects': json.loads((P / 'source-contract.json').read_text())['native_objects'],
    'scope': 'One private provider build only. No loading/execution/FP equivalence/benefit/scientific/release qualification or adoption.',
}, indent=2) + '\n')
print(json.dumps({'provider_build_passed': True, 'library_identity': pin(library)}))
