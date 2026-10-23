"""One incremental provider build; retain every sealed, reused build input."""
from pathlib import Path
import hashlib, json, shutil, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
P = R / 'bin/linux-gcc/lavapipe-intergroup-row-owner/provider'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
plan = doc(W / 'incremental-plan.json')
assert plan['source_independently_reviewed']
assert not P.exists()
for key, value in plan['reused_build_inputs'].items(): assert pin(R / key) == value, key
assert pin(B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so') == plan['previous_target']
subprocess.run(['ninja', '-C', str(B), '-j2',
                'src/gallium/targets/lavapipe/libvulkan_lvp.so'], check=True)
for key, value in plan['reused_build_inputs'].items(): assert pin(R / key) == value, key
target = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'
P.mkdir(parents=True); driver = P / 'libvulkan_lvp.so'; shutil.copy2(target, driver)
assert pin(driver) == pin(target) and pin(driver) != plan['previous_target']
icd = P / 'lvp_icd.json'
icd.write_text(json.dumps({'file_format_version': '1.0.0',
    'ICD': {'library_path': str(driver), 'api_version': '1.4.305'}}, indent=2) + '\n')
result = {'passed': True, 'library': str(driver.relative_to(R)), 'library_identity': pin(driver),
    'ICD': str(icd.relative_to(R)), 'ICD_identity': pin(icd),
    'actual_build_tree_target': str(target.relative_to(R)), 'actual_build_tree_target_identity': pin(target),
    'reused_build_inputs_unchanged': len(plan['reused_build_inputs']),
    'expected_mutable_build_outputs': plan['mutable_outputs'],
    'new_Sirius_CPP_or_shader_compilations': 0, 'provider_builds': 1,
    'global_installation_performed': False, 'provider_or_source_adopted': False,
    'scope': 'One incremental private provider build with unchanged original configuration and sealed reused artifacts. Staged actual provider bytes are separate from mutable compiler workspace. Compilation is not numerical/performance/scientific/release qualification.'}
(W / 'provider-build-result.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'provider_build_passed': True, 'staged_provider': result['library_identity'],
                  'reused_build_inputs_unchanged': result['reused_build_inputs_unchanged']}))
