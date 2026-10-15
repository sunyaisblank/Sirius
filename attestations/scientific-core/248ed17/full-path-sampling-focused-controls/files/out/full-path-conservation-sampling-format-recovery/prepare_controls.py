"""Prepare four original sampling controls only after the current strict producers pass.

Root-operated metadata inspection and immutable input capture; this does not run
tests, compile, validate SPIR-V, or grant scientific/release qualification.
"""
from pathlib import Path
import ast
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import struct
import sys

REVISION = '248ed178cf1b5c445e8b1a3ff2a67b4f907ee6ff'
ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
BUILD = ROOT / 'bin/linux-gcc'
DEST = WORK / 'controls'
ELF = BUILD / 'tests/backend/sirius_backend_tests'
CASES = [
    ('quartic-event-fractions', 'CpuTraceBoundary.RetainedQuarticEventsKeepOriginalFractionsAndExcludeLaterDisk'),
    ('quartic-volume-rollback', 'CpuTraceBoundary.RetainedQuarticVolumeFailurePublishesNothingAndRecovers'),
    ('retained-dopri-sampler', 'RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates'),
    ('near-extremal-conservation', 'FullPathAcceptance.VulkanRetainedNearExtremalRayConservesIndependentInvariants'),
]
SECONDS = {'quartic-event-fractions': 300, 'quartic-volume-rollback': 300,
           'retained-dopri-sampler': 900, 'near-extremal-conservation': 28800}
DEVICE_CASES = {'retained-dopri-sampler', 'near-extremal-conservation'}
READER = WORK / 'originals/verify-producer-9979.py'
TOOLS = ['/usr/bin/python3.12', '/usr/bin/ctest', '/usr/bin/cmake',
         '/usr/bin/ninja', '/usr/bin/readelf', '/usr/bin/ldd', '/usr/bin/git',
         '/usr/bin/g++-14', '/usr/bin/gcc-14', '/usr/bin/as', '/usr/bin/ld']


def query(argv):
    result = subprocess.run(argv, cwd=ROOT, capture_output=True, timeout=30)
    assert result.returncode == 0, (argv, result.returncode, result.stderr[-4096:])
    return result.stdout


def key(path):
    path = Path(path)
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def identity(path):
    path = Path(path)
    resolved = path.resolve(strict=True)
    assert resolved.is_file(), str(path)
    hasher = hashlib.sha256()
    with resolved.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            hasher.update(chunk)
    return {'resolved_path': str(resolved), 'bytes': resolved.stat().st_size,
            'sha256': hasher.hexdigest()}


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def source():
    return {'revision': query(['/usr/bin/git', 'rev-parse', 'HEAD']).decode().strip(),
            'status': query(['/usr/bin/git', 'status', '--porcelain', '--untracked-files=normal']).decode().strip()}


def terminal_owner(case, action):
    folder = WORK / case
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['source'] == {'revision': REVISION, 'status': ''}
    assert owner['case'] == case and owner['argv'] == action['argv']
    assert owner['status'] == 'terminal' and owner['passed'] is True
    assert owner['child_exit'] == 0 and owner['stop_reason'] is None
    assert owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == []
    assert owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True
    assert owner['observed_births'] and all(row['absent'] is True for row in owner['observed_births'])
    assert owner['limits'] == action['limits']
    assert owner['wall_seconds'] < action['limits']['seconds']
    assert owner['peak_sampled_rss_kib'] <= action['limits']['rss_kib']
    assert owner['minimum_sampled_available_kib'] >= action['limits']['reserve_kib']
    before = json.loads((folder / 'inputs-before.json').read_text())
    assert before == json.loads((folder / 'inputs-after.json').read_text())
    for path in [ROOT / 'tests/backend/full_path_acceptance_test.cpp', WORK / 'actions.json',
                 WORK / 'run_owned.py', WORK / 'owned_guard.py']:
        assert before[key(path)] == identity(path)
    return owner


def main():
    assert Path.cwd() == ROOT and not DEST.exists()
    assert source() == {'revision': REVISION, 'status': ''}
    for name, pin in json.loads((WORK / 'postbuild-source-freeze.json').read_text())['pins'].items():
        assert identity(Path(name)) == pin, name
    # Reject inherited environment before any backend discovery or inspection.
    forbidden = [name for name, value in os.environ.items() if value and
                 (name.startswith(('SIRIUS_', 'VK_', 'LP_', 'MESA_', 'GALLIVM_', 'GTEST_', 'CTEST_'))
                  or name in {'LD_PRELOAD', 'LD_AUDIT', 'LD_LIBRARY_PATH'})]
    assert not forbidden, forbidden
    original_actions = json.loads((WORK / 'actions.json').read_text())
    assert list(original_actions) == ['configure', 'build-backend']
    for action in original_actions.values():
        assert set(action) == {'revision', 'argv', 'kind', 'expected_exit', 'limits', 'frozen_inputs', 'scope', 'software_provider'}
        assert type(action['expected_exit']) is int and action['expected_exit'] == 0
    owners = {name: terminal_owner(name, original_actions[name]) for name in original_actions}
    assert '-DSIRIUS_ATTESTATION_ROOT=' in owners['configure']['argv']
    assert owners['build-backend']['argv'] == ['cmake', '--build', '--preset', 'linux-gcc',
        '--parallel', '2', '--target', 'SiriusSourceGovernance', 'sirius_backend_tests']

    cache = {}
    for line in (BUILD / 'CMakeCache.txt').read_text().splitlines():
        match = re.fullmatch(r'([^/#][^:]*):[^=]*=(.*)', line)
        if match:
            cache[match[1]] = match[2]
    strict = {'CMAKE_BUILD_TYPE': 'Release', 'BUILD_TESTS': 'ON',
        'CMAKE_EXPORT_COMPILE_COMMANDS': 'ON', 'SIRIUS_CONTRACT_MODE': '2',
        'SIRIUS_MANDATORY_TESTS': 'ON', 'SIRIUS_WERROR': 'ON',
        'SIRIUS_REQUIRE_VULKAN_RUNTIME': 'ON', 'SIRIUS_BUILD_VIEWER': 'ON',
        'SIRIUS_SANITIZERS': 'none', 'SIRIUS_ALIGNMENT_MODE': 'qualification',
        'SIRIUS_ATTESTATION_ROOT': '', 'CMAKE_CXX_FLAGS': '', 'CMAKE_C_FLAGS': '',
        'CMAKE_CXX_FLAGS_RELEASE': '-O3 -DNDEBUG', 'CMAKE_C_FLAGS_RELEASE': '-O3 -DNDEBUG'}
    assert {name: cache.get(name) for name in strict} == strict
    receipt_path = BUILD / 'generated/sirius/alignment_receipt.json'
    receipt = json.loads(receipt_path.read_text())
    domains = ['macos-moltenvk', 'macos-native-build', 'physical-imax-5616x4096',
        'physical-radeon-780m', 'viewer-native-window-input', 'windows-native-build',
        'windows-native-vulkan', 'wsl2-dozen']
    assert receipt['source_revision'] == REVISION and receipt['source_tree_clean'] is True
    assert receipt['alignment_mode'] == 'qualification' and receipt['admitted_domains'] == []
    assert receipt['records'] == [] and receipt['pending_domains'] == receipt['required_domains'] == domains
    assert receipt['ultimate_ideal'] == {'required': True, 'satisfied': False, 'state': 'attestation_required'}
    assert receipt['operating_model']['sha256'] == identity(ROOT / 'tests/operating_model.json')['sha256']
    configure_stdout = (WORK / 'configure/stdout.log').read_text()
    assert 'Release packaging: DISABLED' in configure_stdout

    commands = json.loads((BUILD / 'compile_commands.json').read_text())
    entries = [row for row in commands if Path(row['file']).resolve() == ROOT / 'tests/backend/full_path_acceptance_test.cpp']
    assert len(entries) == 1
    compilation = entries[0]
    tokens = compilation.get('arguments') or shlex.split(compilation['command'])
    assert Path(tokens[0]).resolve() == Path('/usr/bin/g++-14').resolve()
    required = {'-DSIRIUS_CONTRACT_MODE=2', '-DSIRIUS_HAS_VULKAN_BACKEND=1',
        '-DSIRIUS_TEST_HAS_KERNELS=1', '-DSIRIUS_HAS_RETAINED_COMPUTE=1',
        '-DSIRIUS_TEST_HAS_RETAINED_CAMERA=1', '-DSIRIUS_RELEASE_RESOURCE_LOCKED=1',
        '-Werror', '-Wall', '-Wextra', '-Wpedantic', '-std=c++2c', '-O3', '-DNDEBUG'}
    assert required <= set(tokens)
    assert '-include' + str(ROOT / 'tests/support/fail_on_skip.h') in tokens
    unsafe = {'-Ofast', '-ffast-math', '-funsafe-math-optimizations', '-fassociative-math',
        '-freciprocal-math', '-fno-signed-zeros', '-ffinite-math-only', '-fno-trapping-math'}
    first_party = [row for row in commands if str(Path(row['file']).resolve()).startswith(str(ROOT / 'src/sirius') + '/')]
    for row in [*first_party, compilation]:
        args = row.get('arguments') or shlex.split(row['command'])
        assert not set(args) & unsafe, row['file']
        assert '-DSIRIUS_CONTRACT_MODE=2' in args and '-Werror' in args, row['file']
        assert not any(arg.startswith('-fsanitize=') for arg in args)
    object_path = tokens[tokens.index('-o') + 1]
    assert not Path(object_path).is_absolute()
    dependency_text = query(['/usr/bin/ninja', '-C', str(BUILD), '-t', 'deps', object_path]).decode()
    dep_lines = [line for line in dependency_text.splitlines() if not line.startswith('ninja: Entering directory ')]
    assert dep_lines and dep_lines[0].endswith(' (VALID)'), dep_lines[:1]
    dependencies = [Path(line.strip()) if Path(line.strip()).is_absolute() else BUILD / line.strip()
                    for line in dep_lines[1:] if line.strip()]
    assert ROOT / 'tests/backend/full_path_acceptance_test.cpp' in dependencies
    assert ROOT / 'tests/support/fail_on_skip.h' in dependencies
    assert ROOT / 'tests/support/separated_geodesic_reference.h' in dependencies

    ninja = (BUILD / 'build.ninja').read_text()
    link = next(line for line in ninja.splitlines() if line.startswith('build tests/backend/sirius_backend_tests '))
    assert 'CXX_EXECUTABLE_LINKER__sirius_backend_tests_Release ' in link
    assert object_path in link and 'sirius_backend_tests_input_volume' in link
    link_tokens = shlex.split(link.split(': ', 1)[1].split(' || ', 1)[0])[1:]
    link_files = [Path(item) if item.startswith('/') else BUILD / item for item in link_tokens if item != '|']
    # Bind all23 staged resources using the actual qualification-copy recipe.
    resources = []
    volume_edge = next(i for i, line in enumerate(ninja.splitlines())
                       if line.startswith('build CMakeFiles/sirius_backend_tests_input_volume '))
    volume_command = ninja.splitlines()[volume_edge + 1]
    assert volume_command.startswith('  COMMAND = ')
    volume_tokens = shlex.split(volume_command.split(' = ', 1)[1])
    for i, token in enumerate(volume_tokens):
        if token != 'copy_if_different':
            continue
        end = next((j for j in range(i + 1, len(volume_tokens)) if volume_tokens[j] == '&&'), len(volume_tokens))
        operands = [Path(value) for value in volume_tokens[i + 1:end]]
        assert len(operands) >= 2
        destination = operands[-1]
        assert len(operands) == 2 or destination.is_dir()
        for canonical in operands[:-1]:
            staged = destination / canonical.name if destination.is_dir() else destination
            assert canonical.is_absolute() and staged.is_absolute()
            assert staged.is_relative_to(ELF.parent / 'resources')
            assert identity(canonical)['bytes'] == identity(staged)['bytes']
            assert identity(canonical)['sha256'] == identity(staged)['sha256']
            resources.append({'canonical': key(canonical), 'staged': key(staged), 'identity': identity(staged)})
    assert len(resources) == 23 and len({row['staged'] for row in resources}) == 23
    assert {key(p) for p in (ELF.parent / 'resources').rglob('*') if p.is_file()} == {row['staged'] for row in resources}

    # Only accepted pure definitions; the archived footer/tool calls are excluded.
    reader_tree = ast.parse(READER.read_text())
    definitions = ast.Module(body=[node for node in reader_tree.body
        if isinstance(node, (ast.Import, ast.ImportFrom)) or
        (isinstance(node, ast.FunctionDef) and node.name in {'arrays', 'elf_symbols', 'elf_arrays'})], type_ignores=[])
    reader = {}
    exec(compile(definitions, 'accepted_9979_pure_ELF_definitions', 'exec'), reader)
    retained = BUILD / 'src/sirius/backend/retained'
    expected_arrays = reader['arrays'](retained / 'retained_kernels.h')
    actual_arrays, _ = reader['elf_arrays'](ELF, expected_arrays)
    assert len(actual_arrays) == 40
    assert len([name for name in expected_arrays if name.endswith('Program')]) == 7
    owner_arrays = {}
    for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',
                                       (retained / 'retained_owner_kernels.h').read_text(), re.S):
        words = [int(word) for word in re.findall(r'(\d+)u', body)]
        assert len(words) == int(count)
        owner_arrays[name] = struct.pack('<' + str(len(words)) + 'I', *words)
    assert set(owner_arrays) == {'kTransportShader', 'kEndpointShader'}
    actual_owner_arrays, _ = reader['elf_arrays'](ELF, owner_arrays, 'retained_owner')

    # Actual source/object/Ninja dependencies, not a historical map replay.
    relevant_sources = [
        'tests/backend/full_path_acceptance_test.cpp', 'tests/backend/retained_dopri_test.cpp',
        'tests/backend/cpu_trace_boundary_test.cpp', 'src/sirius/backend/retained_compute.cpp',
        'src/sirius/backend/retained_integrator.cpp', 'src/sirius/backend/retained_trace_executor.cpp',
        'src/sirius/backend/cpu/geodesic_tracer.cpp', 'src/sirius/core/geodesic_integrator.cpp']
    producer_objects = {}
    producer_files = set()
    producer_dependency_text = {}
    for relative in relevant_sources:
        matches = [row for row in commands if Path(row['file']).resolve() == ROOT / relative]
        assert len(matches) == 1, relative
        entry = matches[0]
        args = entry.get('arguments') or shlex.split(entry['command'])
        assert Path(args[0]).resolve() == Path('/usr/bin/g++-14').resolve()
        assert {'-DSIRIUS_CONTRACT_MODE=2', '-Werror', '-O3', '-DNDEBUG'} <= set(args)
        assert not set(args) & unsafe and not any(arg.startswith('-fsanitize=') for arg in args)
        object_name = args[args.index('-o') + 1]
        text = query(['/usr/bin/ninja', '-C', str(BUILD), '-t', 'deps', object_name]).decode()
        lines = [line for line in text.splitlines() if not line.startswith('ninja: Entering directory ')]
        assert lines and lines[0].endswith(' (VALID)'), relative
        deps = [Path(line.strip()) if Path(line.strip()).is_absolute() else BUILD / line.strip()
                for line in lines[1:] if line.strip()]
        assert ROOT / relative in deps
        if relative == 'src/sirius/backend/retained_compute.cpp':
            assert all(retained / name in deps for name in
                       ('retained_kernels.h', 'retained_owner_kernels.h', 'retained_owner_layout.h'))
        if relative.startswith('tests/'):
            assert ROOT / 'tests/support/fail_on_skip.h' in deps
            assert '-include' + str(ROOT / 'tests/support/fail_on_skip.h') in args
            assert object_name in link_tokens
        obj = BUILD / object_name
        ledger = [line.split('\t') for line in (BUILD / '.ninja_log').read_text().splitlines()
                  if '\t' + object_name + '\t' in line]
        assert ledger and int(ledger[-1][2]) == obj.stat().st_mtime_ns
        producer_files.update([obj, *deps])
        producer_dependency_text[relative] = text
        producer_objects[relative] = {'compiler_entry': entry, 'object': identity(obj),
            'source': identity(ROOT / relative), 'dependency_count': len(deps),
            'dependencies': {key(path): identity(path) for path in deps}, 'ninja_log_row': ledger[-1]}
    build_stdout = (WORK / 'build-backend/stdout.log').read_text()
    assert 'Verifying test-label, operating-model, and repository-structure governance' in build_stdout
    assert 'CTest label file is up to date' in build_stdout
    assert 'Building CXX object tests/backend/CMakeFiles/sirius_backend_tests.dir/full_path_acceptance_test.cpp.o' in build_stdout
    assert 'Linking CXX executable tests/backend/sirius_backend_tests' in build_stdout
    with ELF.open('rb') as stream:
        assert stream.read(6) == b'\x7fELF\x02\x01'
    elf_header = query(['/usr/bin/readelf', '-h', str(ELF)]).decode()
    assert 'Advanced Micro Devices X86-64' in elf_header
    icd = Path('/usr/share/vulkan/icd.d/lvp_icd.json')
    provider = Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')
    declared_provider = json.loads(icd.read_text())['ICD']['library_path']
    # This system ICD uses the loader search name; runtime owners separately
    # require the actual absolute LVP mapping, not an inferred loaded provider.
    assert declared_provider in [provider.name, str(provider)]
    dynamic = {}
    libraries = set()
    for target in [ELF, provider]:
        output = query(['/usr/bin/ldd', str(target)]).decode()
        assert 'not found' not in output
        paths = {Path(item) for item in re.findall(r'(?:=>\s+|^\s*)(/\S+)', output, re.M)}
        assert paths, output
        libraries.update(paths)
        dynamic[key(target)] = output
    assert any(path.name.startswith('libvulkan.so') for path in libraries)

    inventory_raw = query(['/usr/bin/ctest', '--test-dir', str(BUILD), '--show-only=json-v1'])
    inventory = json.loads(inventory_raw)
    registration = {}
    for case, name in CASES:
        matches = [row for row in inventory['tests'] if row['name'] == name]
        assert len(matches) == 1, name
        row = matches[0]
        properties = {prop['name']: prop['value'] for prop in row.get('properties', [])}
        assert 'Mandatory' in properties.get('LABELS', []), (name, properties)
        assert row['command'] == [str(ELF), '--gtest_filter=' + name, '--gtest_also_run_disabled_tests']
        # No original per-case TIMEOUT is declared. Outer walls are diagnostics.
        assert 'TIMEOUT' not in properties, (name, properties)
        if case in DEVICE_CASES:
            assert properties['RESOURCE_LOCK'] in ['sirius_vulkan_device', ['sirius_vulkan_device']]
        registration[case] = row
    list_raw = query([str(ELF), '--gtest_list_tests', '--gtest_filter=' + ':'.join(name for _, name in CASES)])
    listed = []
    suite = None
    for line in list_raw.decode().splitlines():
        if line and not line.startswith(' ') and line.rstrip().endswith('.'):
            suite = line.strip()
        elif line.startswith('  ') and suite:
            listed.append(suite + line.strip().split(' #', 1)[0])
    assert sorted(listed) == sorted(name for _, name in CASES)

    files = {Path(__file__).resolve(), WORK / 'actions.json', WORK / 'run_owned.py', WORK / 'owned_guard.py',
        WORK / 'postbuild-source-freeze.json', WORK / 'postbuild-source-contracts.json',
        WORK / 'prepare-controls.patch', WORK / 'originals/prepare-controls-4498.py', READER,
        BUILD / 'CMakeCache.txt', BUILD / 'compile_commands.json',
        BUILD / 'build.ninja', BUILD / 'CMakeFiles/rules.ninja', receipt_path, ELF, icd, provider,
        BUILD / object_path, BUILD / 'tests/evidence_governance.cmake',
        BUILD / 'generated/sirius/base/operating_model_embedded.h', *map(Path, TOOLS),
        *dependencies, *link_files, *producer_files, *libraries, BUILD / '.ninja_log', BUILD / '.ninja_deps',
        *(retained / name for name in ('retained_kernels.h', 'retained_owner_kernels.h', 'retained_owner_layout.h'))}
    files.update(Path(name) for name in json.loads((WORK / 'postbuild-source-freeze.json').read_text())['pins'])
    # Capture the entire actual registration chain, never governance's output inventory or Testing/.
    files.update(BUILD.rglob('CTestTestfile.cmake'))
    files.update((BUILD / 'tests').rglob('*_tests.cmake'))
    files.update((BUILD / 'tests').rglob('*_include.cmake'))
    files.update((BUILD / 'generated/sirius').rglob('*.h'))
    source_files = query(['/usr/bin/git', 'ls-files', '-z', '--', 'CMakeLists.txt', 'CMakePresets.json',
        'cmake', 'src', 'tests', 'scripts', '.github/workflows/ci.yml']).split(b'\0')
    files.update(ROOT / os.fsdecode(path) for path in source_files if path)
    for row in resources:
        files.update([ROOT / row['canonical'], ROOT / row['staged']])
    for phase in owners:
        files.update(WORK / phase / name for name in ['owner.json', 'inputs-before.json', 'inputs-after.json', 'stdout.log', 'stderr.log'])
    for name in ['SIRIUS_SLANGC', 'SIRIUS_SPIRV_VAL', 'SIRIUS_SPIRV_AS', 'SIRIUS_SPIRV_DIS', 'SIRIUS_SPIRV_OPT']:
        files.add(Path(cache[name]))
    files.add(Path('/opt/slang/lib/libslang.so'))
    pins = {key(path): identity(path) for path in sorted(files)}
    assert not any('source-governance-ctest.json' in path or '/Testing/' in path for path in pins)
    assert source() == {'revision': REVISION, 'status': ''}
    DEST.mkdir()
    for name in ['run_owned.py', 'owned_guard.py']:
        shutil.copyfile(WORK / name, DEST / name)
        assert (DEST / name).read_bytes() == (WORK / name).read_bytes()
    (DEST / 'ctest-registration.json').write_bytes(inventory_raw)
    (DEST / 'gtest-list.txt').write_bytes(list_raw)
    (DEST / 'compiler-dependencies.txt').write_text(dependency_text)
    write(DEST / 'producer-dependencies.json', producer_dependency_text)
    write(DEST / 'producer-objects.json', producer_objects)
    write(DEST / 'embedded-retained-arrays.json', {'actual_main_arrays': actual_arrays,
        'actual_owner_arrays': actual_owner_arrays, 'main_count': 40, 'owner_count': 2,
        'program_tables': 7, 'pure_reader_source': identity(READER),
        'scope': 'Exact generated-header bytes in actual rebuilt ELF. No new SPV or LLVM/native graph validation, numerical execution or qualification.'})
    write(DEST / 'dynamic-dependencies.json', dynamic)
    readiness = {'ready_for_root_review': True, 'source': source(), 'strict_cache': strict,
        'alignment_receipt': receipt, 'build_owners': {name: {'wall_seconds': owner['wall_seconds'],
            'peak_sampled_rss_kib': owner['peak_sampled_rss_kib']} for name, owner in owners.items()},
        'compiler_entry': compilation, 'compiler_f_flags': [arg for arg in tokens if arg.startswith('-f')],
        'compiled_header_dependency_count': len(dependencies), 'backend_link_edge': link,
        'resources': resources, 'actual_registration': registration, 'elf_header': elf_header,
        'retained_arrays': {'main_count': 40, 'owner_count': 2, 'program_tables': 7,
            'receipt': identity(DEST / 'embedded-retained-arrays.json')},
        'producer_object_count': len(producer_objects),
        'producer_objects': identity(DEST / 'producer-objects.json'),
        'cache_policy': 'default/unset; no LP/MESA tuning; LVP/device0 only for device cases',
        'claim_limits': ['No control has run during preparation.', '0/8 external domains admitted; release packaging disabled.',
            'Four unchanged original controls; target no declared per-test TIMEOUT and executor maximum_submission_ms remains0.',
            'CPU controls retain supplied independent quartic/refusal contracts; sampler and target require actual device evidence.',
            '40main+2owner arrays/7program tables bound to current headers/ELF; original40physical/4columns/7RHS/3RK/4quartics assertions unchanged.',
            'External ownership walls28800/900/300s are diagnostic only; no fullscience/native/release or repaired702 success before actual execution.']}
    write(DEST / 'readiness.json', readiness)
    commands_out = {'order': [case for case, _ in CASES], 'stop_first_failure': True,
        'commands': [['/usr/bin/python3.12', '-B', str(DEST / 'run_owned.py'), case] for case, _ in CASES],
        'before_first_case_require': ['Root independently matches every source-freeze identity and admits actual ELF/resources/registration/compiler/action closure.',
            'All four output folders absent; copied ordinary owner/guard byte-exact; no other workload.'],
        'before_next_case_require': ['Previous owner terminalPASS/exit0/no guard, cleanup, remaining birth, source or whole-input change.',
            'Exactly selected original case XML tests1/F/E/S/D0/statusrun/resultcompleted; original assertions passed.',
            'Both Vulkan cases require original observed LVP birth/session mapping; default cache remains unset.']}
    write(DEST / 'commands.json', commands_out)
    generated = [path for path in DEST.iterdir() if path.is_file()]
    frozen = sorted([*pins, *(key(path) for path in generated), key(DEST / 'source-freeze.json')])
    actions = {}
    for index, (case, name) in enumerate(CASES):
        previous = [key(DEST / prior / output) for prior, _ in CASES[:index]
                    for output in ['owner.json', 'inputs-before.json', 'inputs-after.json', 'gtest.xml']]
        actions[case] = {'revision': REVISION, 'argv': registration[case]['command'], 'kind': 'test',
            'software_provider': case in DEVICE_CASES, 'disable_software_shader_cache': False,
            'expected_exit': 0, 'limits': {'seconds': SECONDS[case], 'rss_kib': 4194304, 'reserve_kib': 2097152},
            'frozen_inputs': [*frozen, *previous],
            'scope': name + '; unchanged original assertions/fixtures/no added product deadline or submission guard; outer wall is prospective diagnostic only.',
            'root_required_xml': {'suite': name.split('.', 1)[0], 'name': name.split('.', 1)[1],
                'tests': 1, 'failures': 0, 'errors': 0, 'disabled': 0, 'skipped': 0,
                'status': 'run', 'result': 'completed'}}
        assert type(actions[case]['expected_exit']) is int and actions[case]['expected_exit'] == 0
        assert len(actions[case]['frozen_inputs']) == len(set(actions[case]['frozen_inputs']))
    write(DEST / 'actions.json', actions)
    pins.update({key(path): identity(path) for path in [*generated, DEST / 'actions.json']})
    write(DEST / 'source-freeze.json', {'revision': REVISION, 'files': pins,
        'must_match_before_first_owner': True, 'previous_case_outputs_frozen_only_by_later_actions': True})
    print(json.dumps({'ready_for_root_review': True, 'controls': str(DEST), 'static_pins': len(pins),
        'source_freeze': identity(DEST / 'source-freeze.json'), 'actions': identity(DEST / 'actions.json'),
        'runtime_cases_executed': 0}, indent=2))


if __name__ == '__main__':
    main()
