from pathlib import Path
import hashlib, json, shlex, shutil, subprocess

R = Path.cwd()
W = R / 'out/vulkan-descriptor-reuse'
B = R / 'bin/linux-gcc'
T = B / 'vulkan-descriptor-reuse'
assert not T.exists()
T.mkdir()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
baseline = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
limits = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
document = lambda p: json.loads(p.read_text())
write = lambda p, value: p.write_text(json.dumps(value, indent=2) + '\n')
old = R / 'attestations/alignment/4bcd987/gcc-revision-synchronisation/files/out/gcc-revision-synchronisation'
for name in ['run_owned.py', 'owned_guard.py']:
    shutil.copyfile(old / name, W / name)
tools = [Path(x) for x in ['/usr/bin/g++-14', '/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus',
         '/usr/libexec/gcc/x86_64-linux-gnu/14/collect2', '/usr/bin/x86_64-linux-gnu-ld.bfd']]
generated = [B / 'compile_commands.json', B / 'CMakeCache.txt',
             B / 'src/sirius/backend/retained/retained_kernels.h',
             B / 'tests/backend/retained_camera/program_fixture.h']
assert all(p.is_file() for p in tools + generated)
commands = document(B / 'compile_commands.json')
outputs = {}
actions = {}
for name in ['vulkan_device', 'vulkan_smoke_test', 'retained_compute_test', 'retained_dopri_test']:
    command = next(x for x in commands if Path(x['file']).name == name + '.cpp')
    argv = shlex.split(command['command'])
    destination = T / (name + '.o')
    argv[argv.index('-o') + 1] = str(destination)
    argv += ['-MMD', '-MF', str(T / (name + '.d'))]
    outputs[name] = destination
    actions['compile-' + name] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
        'argv': argv, 'limits': limits, 'frozen_inputs': [str(p) for p in tools + generated],
        'scope': 'Isolated strict GCC object using current candidate source/header; primary 4bcd987 products and qualification receipts remain unchanged.'}
line = (W / 'primary-link-commands.txt').read_text().splitlines()[-1]
argv = shlex.split(line.split('&&')[1])
assert argv[0] == '/usr/bin/g++-14'
executable = T / 'original_candidate_controls'
argv[argv.index('-o') + 1] = str(executable)
for i, word in enumerate(argv):
    if word.endswith('.o'):
        name = Path(word).name.removesuffix('.cpp.o')
        argv[i] = str(outputs.get(name, B / word))
    elif word.endswith(('.a', '.so')) and not word.startswith('/'):
        argv[i] = str(B / word)
argv.insert(3, str(outputs['vulkan_device']))
argv += ['-pthread', '-Wl,-Map=' + str(T / 'candidate.map')]
inputs = [Path(x) for x in argv if x.endswith(('.o', '.a', '.so'))]
actions['link-candidate'] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
    'argv': argv, 'limits': limits, 'frozen_inputs': [str(p) for p in tools + inputs],
    'scope': 'Matched fresh concrete Vulkan adapter and all three concrete-header consumers. Reuse only byte-exact generic 4bcd987 objects/libraries; private link does not confer current native or release qualification.'}
shutil.copytree(B / 'tests/backend/resources', T / 'resources')
resource_inputs = sorted((T / 'resources').rglob('*'))
resource_inputs = [p for p in resource_inputs if p.is_file()]
provider = [Path(x) for x in ['/usr/lib/x86_64-linux-gnu/libvulkan.so',
            '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so', '/usr/share/vulkan/icd.d/lvp_icd.json']]
filters = ['VulkanBackend.IdenticalKernelWordsReusePipelineAcrossBuffers',
           'VulkanBackend.IndependentPairCompletesDistinctKernelsAndRejectsSharedBuffers',
           'VulkanBackend.DescriptorReusePreservesOrderedBindingsThroughEviction']
original = document(R / 'attestations/software-vulkan/41fe5f2/portable-denominator-preparation/files/out/portable-denominator-preparation/original-controls/actions.json')
filters += [original['control-' + str(i)]['argv'][1].removeprefix('--gtest_filter=') for i in range(1, 9)]
for i, name in enumerate(filters):
    actions['control-' + str(i + 1)] = {'revision': revision, 'kind': 'test', 'expected_exit': 0,
        'argv': [str(executable), '--gtest_filter=' + name], 'limits': limits,
        'software_provider': True, 'frozen_inputs': [str(p) for p in [executable] + provider + resource_inputs],
        'scope': 'Real software Vulkan binding/stream correctness witness or unchanged original retained control; finite evidence only.'}
write(W / 'actions.json', actions)
write(W / 'preparation.json', {'candidate_revision': revision, 'baseline_revision': baseline,
    'baseline_executable': {'path': str(B / 'tests/backend/sirius_backend_tests'),
        'bytes': (B / 'tests/backend/sirius_backend_tests').stat().st_size,
        'sha256': sha(B / 'tests/backend/sirius_backend_tests')},
    'primary_generated_inputs': {str(p.relative_to(R)): {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in generated},
    'resources': {str(p.relative_to(R)): {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in resource_inputs},
    'compile_actions': [k for k in actions if k.startswith('compile-')],
    'link_inputs': [str(p) for p in inputs], 'filters': filters})
print(json.dumps({'revision': revision, 'actions': len(actions), 'resource_files': len(resource_inputs), 'resource_bytes': sum(p.stat().st_size for p in resource_inputs)}))
