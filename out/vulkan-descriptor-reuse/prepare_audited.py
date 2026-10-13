from pathlib import Path
import hashlib, json, shlex, shutil, subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc'
T = B / 'vulkan-descriptor-reuse'
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert revision == '8a27321c53332fb76273556987040e9ac6ab8fa2'
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
write = lambda p, x: p.write_text(json.dumps(x, indent=2) + '\n')
source = R / 'attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation/files/out/portable-indexed-limb-accumulation/original-controls/candidate/retained_compute_test.cpp'
audited = W / 'audited-observer.cpp'
assert not audited.exists()
shutil.copyfile(source, audited)
assert audited.read_bytes() == source.read_bytes()
assert 'executed_module_' in audited.read_text()
header = T / 'baseline-header/sirius/backend/vulkan/vulkan_device.h'
header.parent.mkdir(parents=True)
baseline = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
old_header = subprocess.check_output(['git', 'show', baseline + ':src/sirius/backend/vulkan/vulkan_device.h'])
header.write_bytes(old_header)
assert hashlib.sha256(old_header).hexdigest() == '604c1ca5ca914d73e0e1c7288f832309e39a07b1ac3006b98812d0fa7e80f3de'
generic = subprocess.check_output(['git', 'show', baseline + ':src/sirius/backend/device.h'])
assert generic == (R / 'src/sirius/backend/device.h').read_bytes()
shutil.copyfile(W / 'actions.json', W / 'actions-before-observer-repair.json')
actions = json.loads((W / 'actions.json').read_text())
command = next(x for x in json.loads((B / 'compile_commands.json').read_text()) if Path(x['file']).name == 'retained_compute_test.cpp')
original = shlex.split(command['command'])
dependencies = shlex.split((T / 'retained_compute_test.d').read_text().replace('\\\n', '').split(':', 1)[1])
tools = ['/usr/bin/g++-14', '/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus',
         '/usr/libexec/gcc/x86_64-linux-gnu/14/collect2', '/usr/bin/x86_64-linux-gnu-ld.bfd']
limits = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
object_paths = {}
for mode in ['baseline', 'candidate']:
    output = T / ('audited-' + mode + '-observer.o')
    assert not output.exists()
    argv = list(original)
    argv[argv.index('-c') + 1] = str(audited)
    argv[argv.index('-o') + 1] = str(output)
    if mode == 'baseline': argv.insert(1, '-I' + str(T / 'baseline-header'))
    argv += ['-MMD', '-MF', str(T / ('audited-' + mode + '-observer.d'))]
    object_paths[mode] = output
    actions['compile-audited-' + mode] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
        'argv': argv, 'limits': limits, 'frozen_inputs': list(dict.fromkeys(tools + dependencies +
            [str(audited), str(source), str(header), str(B / 'compile_commands.json'), str(B / 'CMakeCache.txt')])),
        'scope': 'Compile exact original archived benchmark observer. Baseline uses byte-exact4b concrete header matching old adapter; candidate uses current8a header. Generic ABI and full numerical workload unchanged; no native qualification claim.'}
    adapter = (B / 'src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o'
               if mode == 'baseline' else T / 'vulkan_device.o')
    dopri = (B / 'tests/backend/CMakeFiles/sirius_backend_tests.dir/retained_dopri_test.cpp.o'
             if mode == 'baseline' else T / 'retained_dopri_test.o')
    argv = ['/usr/bin/g++-14', '-O3', '-DNDEBUG', str(output), str(dopri), str(adapter)]
    argv += [str(B / x) for x in ['src/sirius/backend/libsirius_backend_cpu.a', 'lib/libgtest_main.a',
             'src/sirius/backend/libsirius_backend.a', 'src/sirius/render/libsirius_render.a',
             'lib/libgtest.a', 'src/sirius/backend/libsirius_backend_cpu.a', 'src/sirius/backend/libsirius_backend.a']]
    argv += ['/usr/lib/x86_64-linux-gnu/libvulkan.so', '-ldl', str(B / 'src/sirius/core/libsirius_core.a'),
             str(B / 'src/sirius/libsirius_base.a'), '-pthread',
             '-Wl,-Map=' + str(T / ('audited-' + mode + '.map')),
             '-o', str(T / ('audited_' + mode + '_controls'))]
    actions['link-audited-' + mode] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0,
        'argv': argv, 'limits': limits, 'frozen_inputs': list(dict.fromkeys(tools + [x for x in argv if x.endswith(('.o', '.a', '.so'))])),
        'scope': 'Exact original two-TU retained observer executable. Explicit matching adapter and DP object, identical shared stock libraries/resources; real54 module-command joins enabled for both producers.'}
write(W / 'actions.json', actions)
write(W / 'observer-repair.json', {'candidate_revision': revision, 'baseline_revision': baseline,
    'initial_collection_scope': 'Four successfully executed slots have32 digests/51 markers, but lack54 required module-command joins. The original reader refused; no performance verdict and no original20+2 admissibility.',
    'first_comparator_failure': 'missing original static expected-readback-layout.json copy; immutable reference restored byte-exact',
    'second_comparator_failure': 'KeyError executed_module_0_0; canonical benchmark lacked original archived observer instrumentation',
    'source_original': {'path': str(source.relative_to(R)), 'bytes': source.stat().st_size, 'sha256': sha(source)},
    'source_copy': {'path': str(audited.relative_to(R)), 'bytes': audited.stat().st_size, 'sha256': sha(audited)},
    'baseline_concrete_header': {'path': str(header.relative_to(R)), 'bytes': header.stat().st_size, 'sha256': sha(header)},
    'repair': 'Compile the same immutable original audited source for both producers with matching concrete headers, then collect one valid frozen ABBA. No product source, workload, gate or numerical tolerance change; original deficient observations retained.'})
print(json.dumps({'prepared': True, 'original_audited_source_sha256': sha(audited), 'compile_actions': ['compile-audited-baseline', 'compile-audited-candidate']}))
