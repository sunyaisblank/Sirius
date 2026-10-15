"""Prepare six existing controls only after the current strict producers pass.

Root-operated metadata inspection and immutable input capture; this does not run
tests, compile, validate SPIR-V, or grant scientific/release qualification.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys

REVISION = '4498cc41168ad29fc4c111819a4cb7ccaa54a8bf'
ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
BUILD = ROOT / 'bin/linux-gcc'
DEST = WORK / 'controls'
ELF = BUILD / 'tests/backend/sirius_backend_tests'
GOVERNANCE = 'OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims'
CASES = [
    ('source-governance', GOVERNANCE),
    ('module-capabilities', 'KernelParity.PrecisionProbeArtifactsCarryOnlyTheirDeclaredFloat64Capability'),
    ('fp32', 'KernelParity.RejectedNeighbourRetriesExactMinimumFp32'),
    ('fp32comp', 'KernelParity.RejectedNeighbourRetriesExactMinimumFp32Comp'),
    ('fp64', 'KernelParity.RejectedNeighbourRetriesExactMinimumFp64'),
    ('ellis-physical-ray', 'KernelParity.EllisTwoSheetTraceCrossesThroatAndMapsTheOppositeSky'),
]
LIMITS = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
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
    for path in [ROOT / 'tests/backend/kernel_parity_test.cpp', WORK / 'actions.json',
                 WORK / 'run_owned.py', WORK / 'owned_guard.py']:
        assert before[key(path)] == identity(path)
    return owner


def main():
    assert Path.cwd() == ROOT and not DEST.exists()
    assert source() == {'revision': REVISION, 'status': ''}
    # Reject inherited environment before any backend discovery or inspection.
    forbidden = [name for name, value in os.environ.items() if value and
                 (name.startswith(('SIRIUS_', 'VK_', 'LP_', 'MESA_', 'GALLIVM_', 'GTEST_', 'CTEST_'))
                  or name in {'LD_PRELOAD', 'LD_AUDIT', 'LD_LIBRARY_PATH'})]
    assert not forbidden, forbidden
    original_actions = json.loads((WORK / 'actions.json').read_text())
    assert list(original_actions) == ['configure', 'build-backend']
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
    entries = [row for row in commands if Path(row['file']).resolve() == ROOT / 'tests/backend/kernel_parity_test.cpp']
    assert len(entries) == 1
    compilation = entries[0]
    tokens = compilation.get('arguments') or shlex.split(compilation['command'])
    assert Path(tokens[0]).resolve() == Path('/usr/bin/g++-14').resolve()
    required = {'-DSIRIUS_CONTRACT_MODE=2', '-DSIRIUS_HAS_VULKAN_BACKEND=1',
        '-DSIRIUS_TEST_HAS_KERNELS=1', '-DSIRIUS_RELEASE_RESOURCE_LOCKED=1',
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
    assert ROOT / 'tests/backend/kernel_parity_test.cpp' in dependencies
    assert ROOT / 'tests/support/fail_on_skip.h' in dependencies
    assert ROOT / 'tests/support/test_resource.h' in dependencies

    ninja = (BUILD / 'build.ninja').read_text()
    link = next(line for line in ninja.splitlines() if line.startswith('build tests/backend/sirius_backend_tests '))
    assert 'CXX_EXECUTABLE_LINKER__sirius_backend_tests_Release ' in link
    assert object_path in link and 'sirius_backend_tests_input_volume' in link
    link_tokens = shlex.split(link.split(': ', 1)[1].split(' || ', 1)[0])[1:]
    link_files = [Path(item) if item.startswith('/') else BUILD / item for item in link_tokens if item != '|']
    shader_dependencies = set()
    resources = []
    for name in ['parity_probe.spv', 'parity_probe_fp32comp.spv', 'parity_probe_fp64.spv']:
        edge = next(line for line in ninja.splitlines() if line.startswith('build kernels/' + name + ' '))
        assert ': CUSTOM_COMMAND ' in edge
        deps = [Path(item) for item in shlex.split(edge.split(': CUSTOM_COMMAND ', 1)[1])]
        assert ROOT / 'src/sirius/kernels/parity_probe.slang' in deps
        assert ROOT / 'src/sirius/kernels/gr_integrator.slang' in deps
        assert ROOT / 'src/sirius/kernels/gr_ellis_topology.slang' in deps
        shader_dependencies.update(deps)
        canonical = BUILD / 'kernels' / name
        staged = ELF.parent / 'resources/kernels' / name
        canonical_identity = identity(canonical)
        staged_identity = identity(staged)
        assert (canonical_identity['bytes'], canonical_identity['sha256']) == (staged_identity['bytes'], staged_identity['sha256'])
        assert f'copy_if_different {canonical} {staged}' in ninja
        resources.append({'canonical': key(canonical), 'staged': key(staged), 'identity': staged_identity})
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
        if case == 'source-governance':
            assert properties['TIMEOUT'] == 30.0
            assert row['command'][0] == '/usr/bin/cmake'
            assert row['command'][-2:] == ['-P', str(BUILD / 'tests/evidence_governance.cmake')]
        else:
            assert row['command'] == [str(ELF), '--gtest_filter=' + name, '--gtest_also_run_disabled_tests']
            lock = properties['RESOURCE_LOCK']
            assert lock in ['sirius_vulkan_device', ['sirius_vulkan_device']]
        registration[case] = row
    list_raw = query([str(ELF), '--gtest_list_tests', '--gtest_filter=' + ':'.join(name for _, name in CASES[1:])])
    listed = []
    suite = None
    for line in list_raw.decode().splitlines():
        if line == 'KernelParity.':
            suite = line
        elif line.startswith('  ') and suite:
            listed.append(suite + line.strip().split(' #', 1)[0])
    assert sorted(listed) == sorted(name for _, name in CASES[1:])

    files = {Path(__file__).resolve(), WORK / 'actions.json', WORK / 'run_owned.py', WORK / 'owned_guard.py',
        WORK / 'root-source-admission.json', BUILD / 'CMakeCache.txt', BUILD / 'compile_commands.json',
        BUILD / 'build.ninja', BUILD / 'CMakeFiles/rules.ninja', receipt_path, ELF, icd, provider,
        BUILD / object_path, BUILD / 'tests/evidence_governance.cmake',
        BUILD / 'generated/sirius/base/operating_model_embedded.h', *map(Path, TOOLS),
        *dependencies, *link_files, *shader_dependencies, *libraries}
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
    write(DEST / 'dynamic-dependencies.json', dynamic)
    readiness = {'ready_for_root_review': True, 'source': source(), 'strict_cache': strict,
        'alignment_receipt': receipt, 'build_owners': {name: {'wall_seconds': owner['wall_seconds'],
            'peak_sampled_rss_kib': owner['peak_sampled_rss_kib']} for name, owner in owners.items()},
        'compiler_entry': compilation, 'compiler_f_flags': [arg for arg in tokens if arg.startswith('-f')],
        'compiled_header_dependency_count': len(dependencies), 'backend_link_edge': link,
        'resources': resources, 'actual_registration': registration, 'elf_header': elf_header,
        'cache_policy': 'default/unset; no LP/MESA tuning; LVP/device0 only for device cases',
        'claim_limits': ['No control has run during preparation.', '0/8 external domains admitted; release packaging disabled.',
            'Three rung controls cover the real adaptive candidate and central controller/state, not all continuation fields or frame publication.',
            'The original Ellis finite physical ray is a separate regression; no Kerr frame/CI-causality or full P1/P6 qualification claim.']}
    write(DEST / 'readiness.json', readiness)
    commands_out = {'order': [case for case, _ in CASES], 'stop_first_failure': True,
        'commands': [['/usr/bin/python3.12', '-B', str(DEST / 'run_owned.py'), case] for case, _ in CASES],
        'before_first_case_require': ['Root independently matches all source-freeze file identities to current bytes and admits readiness/actions.',
            'All six case output folders are absent; exact helpers equal the accepted root originals.'],
        'before_next_case_require': ['owner passed/terminal/exit0/no stop, cleanup or live birth; source and whole inputs unchanged',
            'For kind:test: exactly the selected gtest in XML; tests1/errors0/failures0/disabled0/skipped0, status=run/result=completed.',
            'For each rung: its exact *_evidence property equals real_adaptive_minimum_retry_controls_executed; FP64 unsupported_kernel_declined is refusal-only, not executed FP64 evidence.',
            'For governance: original CTest JUnit contains exactly the selected test with no failure/error/skipped.']}
    write(DEST / 'commands.json', commands_out)
    generated = [path for path in DEST.iterdir() if path.is_file()]
    frozen = sorted([*pins, *(key(path) for path in generated), key(DEST / 'source-freeze.json')])
    actions = {}
    for index, (case, name) in enumerate(CASES):
        argv = (['/usr/bin/ctest', '--test-dir', str(BUILD), '-V', '--output-on-failure', '--no-tests=error',
            '--output-junit', '{folder}/ctest.xml', '--output-log', '{folder}/ctest.log', '-R', '^' + re.escape(name) + '$']
            if case == 'source-governance' else [str(ELF), '--gtest_filter=' + name])
        previous = [key(DEST / prior / output) for prior, _ in CASES[:index]
                    for output in ['owner.json', 'inputs-before.json', 'inputs-after.json',
                        'ctest.xml' if prior == 'source-governance' else 'gtest.xml']]
        actions[case] = {'revision': REVISION, 'argv': argv, 'kind': 'build' if index == 0 else 'test',
            'software_provider': case in {'fp32', 'fp32comp', 'fp64', 'ellis-physical-ray'},
            'disable_software_shader_cache': False, 'expected_exit': 0, 'limits': LIMITS,
            'frozen_inputs': [*frozen, *previous],
            'scope': name + '; unchanged original assertions/fixtures and applicable RunProbe1000ms; finite focused control only.'}
        if index:
            actions[case]['root_required_xml'] = {'suite': 'KernelParity', 'name': name.split('.', 1)[1],
                'tests': 1, 'failures': 0, 'errors': 0, 'disabled': 0, 'skipped': 0,
                'status': 'run', 'result': 'completed'}
        if case in {'fp32', 'fp32comp', 'fp64'}:
            actions[case]['root_required_xml']['property'] = {
                'name': case + '_evidence', 'value': 'real_adaptive_minimum_retry_controls_executed'}
    write(DEST / 'actions.json', actions)
    pins.update({key(path): identity(path) for path in [*generated, DEST / 'actions.json']})
    write(DEST / 'source-freeze.json', {'revision': REVISION, 'files': pins,
        'must_match_before_first_owner': True, 'previous_case_outputs_frozen_only_by_later_actions': True})
    print(json.dumps({'ready_for_root_review': True, 'controls': str(DEST), 'static_pins': len(pins),
        'source_freeze': identity(DEST / 'source-freeze.json'), 'actions': identity(DEST / 'actions.json'),
        'runtime_cases_executed': 0}, indent=2))


if __name__ == '__main__':
    main()
