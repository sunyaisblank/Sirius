"""Relink the original render test and runbook against the qualified backend."""
from pathlib import Path
import hashlib, json, os, re, shlex, shutil, subprocess, sys
R = Path.cwd(); W = Path(__file__).resolve().parent; P = W.parent
C = R / 'bin/linux-gcc'; B = C / 'native-original-input-ownership'
FRAME = B / 'frame'; SOURCE = '08add7e5633bfbb5aa19343e45497cd33a2edaa3'
def doc(p): return json.loads(p.read_text())
def write(p, v): p.write_text(json.dumps(v, indent=2) + '\n')
def pin(p):
    return {'bytes': p.stat().st_size,
            'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
accepted = doc(P / 'readers/controls-1-8-acceptance.json')
assert accepted['all8_original_controls_pass'] and accepted['all22_original_benefit_checks_pass']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == SOURCE
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
mode, = sys.argv[1:]

if mode == 'build-actions':
    assert not FRAME.exists(); FRAME.mkdir()
    target = 'src/sirius/app/sirius_render_tests'
    lines = subprocess.check_output(['ninja', '-C', str(C), '-t', 'commands', target], text=True).splitlines()
    line, = [x for x in lines if ' -o ' + target + ' ' in x]
    args = shlex.split(line); start = args.index('/usr/bin/g++-14'); end = args.index('&&', start)
    args = args[start:end]
    original_inputs = []
    for word in args:
        p = Path(word); p = p if p.is_absolute() else C / p
        if p.is_file(): original_inputs.append(p)
    original_exe = C / target
    assert pin(original_exe) == {'bytes': 12227616,
        'sha256': '662fdeda2b170f21846bccf589f14836dfb69e5df74f30dad377fbc3acd2abdd'}
    assert args.count('src/sirius/backend/libsirius_backend.a') == 1
    args = [str(B / 'matched/libsirius_backend.a') if x == 'src/sirius/backend/libsirius_backend.a' else x for x in args]
    args[args.index('-o') + 1] = str(FRAME / 'sirius_render_tests')
    volume = C / 'src/sirius/app/resources'
    resources = {str(p.relative_to(volume)): pin(p) for p in volume.rglob('*') if p.is_file()}
    assert resources and all(p.resolve().is_relative_to(volume.resolve()) for p in volume.rglob('*') if p.is_file())
    recipe = {'argv': args, 'cwd': str(C), 'original_executable': str(original_exe),
        'original_executable_identity': pin(original_exe), 'original_resources': str(volume),
        'resources': resources, 'backend': pin(B / 'matched/libsirius_backend.a'),
        'source': SOURCE, 'new_CPP_compilations': 0, 'changed_link_archive': 'sirius_backend'}
    write(W / 'link-recipe.json', recipe)
    inputs = set(original_inputs) | {B / 'matched/libsirius_backend.a', Path(__file__).resolve(),
        W / 'link-recipe.json', P / 'readers/controls-1-8-acceptance.json',
        P / 'matched-result.json', P / 'compiler-result.json'}
    inputs |= {p for p in volume.rglob('*') if p.is_file()}
    action = {'revision': SOURCE, 'kind': 'diagnostic', 'expected_exit': 0,
        'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'argv': ['/usr/bin/python3', '-B', str(Path(__file__).resolve()), 'build'],
        'frozen_inputs': sorted(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()) for p in inputs),
        'software_provider': False,
        'scope': 'One link of the original render test with the already-qualified private backend archive; reuse every original render/app/test object and exact original executable-adjacent resources. No C++ compilation, frame execution, source/gate/scene change or release admission.'}
    write(W / 'actions.json', {'build': action})
elif mode == 'build':
    recipe = doc(W / 'link-recipe.json')
    assert pin(B / 'matched/libsirius_backend.a') == recipe['backend']
    subprocess.run(recipe['argv'], cwd=recipe['cwd'], check=True)
    volume = Path(recipe['original_resources']); destination = FRAME / 'resources'
    assert not destination.exists(); shutil.copytree(volume, destination)
    actual = {str(p.relative_to(destination)): pin(p) for p in destination.rglob('*') if p.is_file()}
    assert actual == recipe['resources']
    assert pin(Path(recipe['original_executable'])) == recipe['original_executable_identity']
    write(W / 'build-result.json', {'passed': True, 'source': SOURCE,
        'executable': str((FRAME / 'sirius_render_tests').relative_to(R)),
        'executable_identity': pin(FRAME / 'sirius_render_tests'), 'original_objects_reused': True,
        'original_resources_copied_exactly': actual, 'new_CPP_compilations': 0,
        'original_render_binary_unchanged': True, 'release_readiness_claimed': False})
elif mode == 'frame-actions':
    assert doc(W / 'build/owner.json')['passed'] and doc(W / 'build-result.json')['passed']
    assert not (W / 'frame').exists()
    exe = FRAME / 'sirius_render_tests'
    dependencies = set()
    for obj in [exe, Path('/usr/bin/cmake'), R / doc(P / 'providers.json')['providers']['candidate']['library']]:
        output = subprocess.check_output(['ldd', str(obj)], text=True)
        dependencies |= {Path(x).resolve() for x in re.findall(r'(/[^\s()]+)', output) if Path(x).is_file()}
    # Require the same actual global dependency bytes as the new C8 cohort.
    c8 = doc(P / 'runtime/control-8/inputs-before.json')
    for dependency in dependencies:
        key = str(dependency)
        if key in c8:
            assert pin(dependency) == {k: c8[key][k] for k in ['bytes', 'sha256']}
    inputs = dependencies | {exe,
        P / 'readers/controls-1-8-acceptance.json', P / 'readers/benefit.json',
        P / 'readers/decimal-agreement.json', P / 'runtime/control-8/inputs-before.json',
        P / 'providers.json', P / 'compiler-result.json', P / 'matched-result.json',
        W / 'build-result.json', R / 'tests/operational/vulkan_runtime_profile.cmake',
        Path('/usr/bin/cmake'), Path('/usr/bin/python3'), Path('/etc/ld.so.cache')}
    inputs |= {p for p in W.iterdir() if p.is_file()}
    inputs |= {p for p in FRAME.rglob('*') if p.is_file()}
    provider = doc(P / 'providers.json')['providers']['candidate']
    inputs |= {Path(provider['library']), Path(provider['ICD'])}
    temp = R / 'renders/native-original-input-08add7e5/cold-frame-tmp'
    assert not temp.exists(); temp.mkdir(parents=True)
    actions = doc(W / 'actions.json')
    actions['frame'] = {'revision': SOURCE, 'kind': 'test', 'expected_exit': 0,
        'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'argv': ['/usr/bin/cmake', '-DSIRIUS_RENDER_TESTS=' + str(exe), '-P',
                 str(R / 'tests/operational/vulkan_runtime_profile.cmake')],
        'frozen_inputs': sorted(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()) for p in inputs),
        'software_provider': True, 'runtime_provider': provider['library'], 'runtime_icd': provider['ICD'],
        'build_environment': {'TMPDIR': str(temp)},
        'scope': 'One original cold OperationalVulkan.RequiredRuntimeDispatch/RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent,128x128 one-SPP Kerr with unchanged defaults, cache-disable policy, governors and original output assertions. All8 controls/all22 gates passed first. Original test removes its own temporary PPM after checking it; TMPDIR only places that original disposable output under renders/. All existing global dependency bytes shared with the new C8 cohort are required exact. No altered feature/property/scene/count/budget, retry or full scientific/Mandatory/release claim.'}
    write(W / 'actions.json', actions)
else:
    raise ValueError(mode)
