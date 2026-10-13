"""Prepare one private, matched original-control producer from preserved shaders."""
from pathlib import Path
import hashlib, json, re, shlex, shutil, struct, subprocess

R = Path.cwd().resolve()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc'
T = B / 'guarded-product-high-controls'
REVISION = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == REVISION
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
assert not T.exists()
T.mkdir()
def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def write(name, value):
    (W / name).write_text(json.dumps(value, indent=2) + '\n')
numeric = R / 'attestations/software-vulkan/4bcd987/guarded-product-high-numerics/files/out/guarded-product-high-numerics'
for name in ['run_owned.py', 'owned_guard.py']:
    shutil.copyfile(numeric / name, W / name)
protected_path = R / 'attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation/files/out/portable-indexed-limb-accumulation/protected-before.json'
protected = json.loads(protected_path.read_text())
for name, pin in protected.items():
    assert identity(R / name) == pin, name
write('protected-before.json', protected)
old97 = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse'
oldwork = old97 / 'files/out/vulkan-descriptor-reuse'
prior_link = json.loads((oldwork / 'link-audited-baseline/inputs-before.json').read_text())
reused = {}
for name, pin in prior_link.items():
    if name.endswith(('.o', '.a', '.so')) and 'vulkan-descriptor-reuse' not in name:
        path = R / name
        assert identity(path) == {k: pin[k] for k in ['bytes', 'sha256']}, name
        reused[name] = pin
observer_origin = R / 'attestations/software-vulkan/6014f91/portable-indexed-limb-accumulation/files/out/portable-indexed-limb-accumulation/original-controls/candidate/retained_compute_test.cpp'
observer = W / 'audited-observer.cpp'
assert identity(observer_origin) == {'bytes': 335033, 'sha256': '399cd52cc6d9fc89e4b97ff63b06a54b06b4f80128f5a99c2c6d1278c62f3dc3'}
assert observer_origin.read_bytes() == (oldwork / 'audited-observer.cpp').read_bytes()
shutil.copyfile(observer_origin, observer)
baseline_origin = old97 / 'input-versions/41eb7ce5e62380913a6c70fa2c3228acf8c427831c0611c84a670e238dc32856'
assert identity(baseline_origin) == {'bytes': 10753968, 'sha256': '41eb7ce5e62380913a6c70fa2c3228acf8c427831c0611c84a670e238dc32856'}
baseline = T / 'original_baseline_controls'
shutil.copyfile(baseline_origin, baseline)
baseline.chmod(0o755)
pilot = R / 'attestations/software-vulkan/4bcd987/guarded-product-high-screen'
manifest = json.loads((pilot / 'manifest.json').read_text())
header_origin = B / 'src/sirius/backend/retained/retained_kernels.h'
header = T / 'headers/retained_kernels.h'
header.parent.mkdir()
text = header_origin.read_text()
pattern = re.compile(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', re.S)
def arrays(content):
    records = {}
    for count, name, body in pattern.findall(content):
        words = [int(v) for v in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count), name
        records[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    assert len(records) == 40
    return records
stock_arrays = arrays(text)
replacements = {'kEndpointPortableShader': 'endpoint-core',
                'kEndpointPortableNormalSumShader': 'endpoint-normal',
                'kTransportPortableNormalSumShader': 'transport-normal',
                'kDopriPhasePortableNormalSumShader': 'dopri_phase-normal',
                'kDensePortableShader': 'dense-core'}
replacement_bytes = {}
origins = [protected_path, observer_origin, baseline_origin, pilot / 'manifest.json', oldwork / 'link-audited-baseline/inputs-before.json']
for name, variant in replacements.items():
    for mode in ['stock', 'candidate']:
        relative = 'files/out/guarded-product-high-screen/modules/' + mode + '-' + variant + '.spv'
        path = pilot / relative
        assert identity(path) == {k: manifest['files'][relative][k] for k in ['bytes', 'sha256']}
        origins.append(path)
        if mode == 'stock':
            assert stock_arrays[name] == path.read_bytes(), name
        else:
            replacement_bytes[name] = path.read_bytes()
def replace(match):
    name = match[2]
    if name not in replacement_bytes:
        return match[0]
    data = replacement_bytes[name]
    assert len(data) % 4 == 0
    words = struct.unpack('<' + str(len(data) // 4) + 'I', data)
    lines = [','.join(str(v) + 'u' for v in words[i:i+16]) + ',' for i in range(0, len(words), 16)]
    return 'inline constexpr std::array<std::uint32_t, ' + str(len(words)) + '> ' + name + '{{\n' + '\n'.join(lines) + '\n}};'
header.write_text(pattern.sub(replace, text))
candidate_arrays = arrays(header.read_text())
for name, data in stock_arrays.items():
    assert candidate_arrays[name] == replacement_bytes.get(name, data), name
assert sum(candidate_arrays[name] != data for name, data in stock_arrays.items()) == 5
shutil.copytree(B / 'tests/backend/resources', T / 'resources')
tools = [Path(p) for p in ['/usr/bin/g++-14', '/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus',
         '/usr/libexec/gcc/x86_64-linux-gnu/14/collect2', '/usr/bin/x86_64-linux-gnu-ld.bfd', '/usr/bin/python3',
         '/usr/bin/readelf', '/usr/bin/nm', '/usr/bin/objdump']]
generated = [B / 'compile_commands.json', B / 'CMakeCache.txt', header_origin,
             B / 'src/sirius/backend/retained/retained_owner_kernels.h',
             B / 'src/sirius/backend/retained/retained_owner_layout.h',
             B / 'tests/backend/retained_camera/program_fixture.h', header]
google_headers = [p for p in (B / '_deps/googletest-src/googletest').rglob('*') if p.is_file() and p.suffix in ['.h', '.hpp']]
shared = [*tools, *generated, *google_headers, *origins, W / 'prepare.py', W / 'audited-observer.cpp']
shared = list(dict.fromkeys(str(p.resolve()) for p in shared))
commands = json.loads((B / 'compile_commands.json').read_text())
limits = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
actions = {}
objects = {}
for name, source in [('compute', R / 'src/sirius/backend/retained_compute.cpp'), ('observer', observer)]:
    original_name = 'retained_compute.cpp' if name == 'compute' else 'retained_compute_test.cpp'
    command = next(c for c in commands if Path(c['file']).name == original_name)
    argv = shlex.split(command['command'])
    argv.insert(1, '-I' + str(header.parent))
    argv[argv.index('-c') + 1] = str(source)
    obj = T / (name + '.o')
    argv[argv.index('-o') + 1] = str(obj)
    argv.extend(['-MMD', '-MF', str(T / (name + '.d'))])
    objects[name] = obj
    actions['compile-' + name] = {'revision': REVISION, 'kind': 'diagnostic', 'expected_exit': 0,
        'argv': argv, 'limits': limits, 'frozen_inputs': shared,
        'scope': 'Current strict GCC preset command, private five-array generated header shared by implementation and immutable original audited observer. No primary build writes or full qualification verdict.'}
original_link = json.loads((oldwork / 'actions.json').read_text())['link-audited-baseline']['argv']
argv = list(original_link)
argv[3] = str(objects['observer'])
argv.insert(4, str(objects['compute']))
argv[argv.index('-o') + 1] = str(T / 'original_candidate_controls')
argv = ['-Wl,-Map=' + str(T / 'candidate.map') if s.startswith('-Wl,-Map=') else s for s in argv]
link_inputs = [s for s in argv if s.endswith(('.o', '.a', '.so'))]
actions['link-candidate'] = {'revision': REVISION, 'kind': 'diagnostic', 'expected_exit': 0,
    'argv': argv, 'limits': limits, 'frozen_inputs': list(dict.fromkeys(shared + link_inputs)),
    'scope': 'Explicit matched observer and retained implementation override old archive compute member; byte-exact current stock adapter, DP object and shared libraries. Original forty arrays checked in actual linked ELF before GPU controls.'}
write('actions.json', actions)
write('preparation.json', {'source_revision': REVISION, 'replacements': replacements,
    'stock_header': {'path': str(header_origin.relative_to(R)), **identity(header_origin)},
    'candidate_header': {'path': str(header.relative_to(R)), **identity(header)},
    'array_count': 40, 'changed_arrays': 5, 'unchanged_arrays': 35,
    'observer_origin': {'path': str(observer_origin.relative_to(R)), **identity(observer_origin)},
    'baseline_origin': {'path': str(baseline_origin.relative_to(R)), **identity(baseline_origin)},
    'baseline_copy': {'path': str(baseline.relative_to(R)), **identity(baseline)},
    'current_byte_exact_reused_objects_and_libraries': reused,
    'shared_nontracked_inputs': shared,
    'scope': 'Private original integrated correctness producer preparation; no benefit or adopted production policy.'})
print(json.dumps({'prepared': True, 'arrays': 40, 'changed': 5, 'shared_inputs': len(shared), 'compile_actions': ['compile-compute', 'compile-observer']}))
