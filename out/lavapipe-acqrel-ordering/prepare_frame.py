"""Reuse the exact original render ELF/resources; only the candidate driver differs."""
from pathlib import Path
import hashlib, json, re, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent
OLD = R / 'out/native-original-input-ownership'
F = W / 'cold-frame'; assert not F.exists()
SOURCE = '1c95c98f61696541f520d871f7b492035074089b'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
accepted = doc(W / 'readers/controls-1-8-acceptance.json')
assert accepted['all8_original_controls_pass'] and accepted['all22_original_benefit_checks_pass']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
old_build = doc(OLD / 'cold-frame/build-result.json')
exe = R / old_build['executable']
assert old_build['passed'] and pin(exe) == old_build['executable_identity']
resources = {str(p.relative_to(exe.parent / 'resources')): pin(p)
    for p in (exe.parent / 'resources').rglob('*') if p.is_file()}
assert resources == old_build['original_resources_copied_exactly']
provider = doc(W / 'providers.json')['providers']['candidate']
assert pin(Path(provider['library'])) == provider['library_identity']
assert pin(Path(provider['ICD'])) == provider['ICD_identity']
dependencies = set()
for obj in [exe, Path('/usr/bin/cmake'), Path(provider['library'])]:
    output = subprocess.check_output(['ldd', str(obj)], text=True)
    dependencies |= {Path(x).resolve() for x in re.findall(r'(/[^\s()]+)', output) if Path(x).is_file()}
c8 = doc(W / 'runtime/control-8/inputs-before.json')
for dependency in dependencies:
    if str(dependency) in c8:
        assert pin(dependency) == {k: c8[str(dependency)][k] for k in ['bytes', 'sha256']}
F.mkdir()
(F / 'owned_guard.py').symlink_to('../../native-cpu-row-ownership/cold-frame/owned_guard.py')
(F / 'run_owned.py').write_text('"""Reuse the preserved original frame owner."""\nfrom pathlib import Path\n_preserved_runner = Path(__file__).resolve().parents[2] / "native-cpu-row-ownership/cold-frame/run_owned.py"\nexec(compile(_preserved_runner.read_text(), str(_preserved_runner), "exec"))\n')
reuse = {'passed': True, 'original_build_result': str((OLD / 'cold-frame/build-result.json').relative_to(R)),
    'original_build_result_identity': pin(OLD / 'cold-frame/build-result.json'),
    'executable': str(exe.relative_to(R)), 'executable_identity': pin(exe),
    'original_resources_copied_exactly': resources, 'new_CPP_compilations': 0, 'new_links': 0,
    'scope': 'Prior exact original render executable/resources reused; its prior frame failed. Only the new provider changes this frame attempt.'}
(F / 'reuse.json').write_text(json.dumps(reuse, indent=2) + '\n')
inputs = dependencies | {exe, Path(provider['library']), Path(provider['ICD']), Path(__file__).resolve(),
    R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
    R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
    OLD / 'cold-frame/build-result.json', W / 'providers.json', W / 'provider-identity.json',
    W / 'readers/controls-1-8-acceptance.json', W / 'readers/benefit.json',
    W / 'readers/decimal-agreement.json', W / 'runtime/control-8/inputs-before.json',
    W / 'ordinary/ordinary-code-review.json', R / 'tests/operational/vulkan_runtime_profile.cmake',
    Path('/usr/bin/cmake'), Path('/usr/bin/python3'), Path('/etc/ld.so.cache')}
inputs |= {p for p in F.iterdir() if p.is_file()}
inputs |= {p for p in (exe.parent / 'resources').rglob('*') if p.is_file()}
temp = R / 'renders/lavapipe-acqrel-1c95c98f/cold-frame-tmp'
assert not temp.exists(); temp.mkdir(parents=True)
action = {'revision': SOURCE, 'kind': 'test', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'argv': ['/usr/bin/cmake', '-DSIRIUS_RENDER_TESTS=' + str(exe), '-P',
             str(R / 'tests/operational/vulkan_runtime_profile.cmake')],
    'frozen_inputs': sorted(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()) for p in inputs),
    'software_provider': True, 'runtime_provider': provider['library'], 'runtime_icd': provider['ICD'],
    'build_environment': {'TMPDIR': str(temp)},
    'scope': 'One changed original cold OperationalVulkan.RequiredRuntimeDispatch/RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent,128x128 one-SPP Kerr. All8/all22 passed first. Same exact original render ELF/resources,0compiles/links; only candidate provider lowering changes. Original scene/defaults/counts/fields/cache-disable/governors/limits/assertions unchanged. Original test removes its checked temporary PPM; TMPDIR only places that original disposable output under renders/. No unchanged retry, altered admission or scientific/Mandatory/release claim.'}
(F / 'actions.json').write_text(json.dumps({'frame': action}, indent=2) + '\n')
print(json.dumps({'exact_existing_render_ELF_reused': True, 'new_compilations': 0, 'new_links': 0,
                  'frame_executed': False, 'source': SOURCE}))
