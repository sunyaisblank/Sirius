"""Replay governed Sirius science cohorts with the ordinary Linux GCC products."""

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
BUILD = ROOT / 'bin/linux-gcc'
WORK = ROOT / 'out/current-integrated-science'
RENDERS = ROOT / 'renders/current-integrated-science'
MODEL = ROOT / 'tests/operating_model.json'
GUARD = ROOT / 'out/native-cpu-row-ownership/cold-frame'
HISTORICAL = ROOT / 'attestations/software-vulkan/b34d8c2/cpu-dense-growth-only-integration/diagnostic'
E2_ACTIONS = ('34a94cd5:attestations/scientific-core/f22dc025/e2-original-eleven-cpu-witnesses/'
              'files/out/e2-original-cpu-witnesses-pending/actions.json')
P4_ACTIONS = ('e86b5baf:attestations/scientific-core/4bcd987/p4-page-thorne-current-source/'
              'files/out/science-p4-page-thorne/actions.json')
ICD = Path('/usr/share/vulkan/icd.d/lvp_icd.json')
PROVIDER = Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')
CTEST = ['/usr/bin/ctest', '--test-dir', str(BUILD), '-C', 'Release']
CHECK_KEYS = ('failures', 'errors', 'skipped', 'disabled')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def pin(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path),
            'bytes': path.stat().st_size, 'sha256': digest}


def read(path):
    return json.loads(path.read_text())


def write(path, data):
    path.write_text(json.dumps(data, sort_keys=True, indent=2) + '\n')


def git_object(revision_path):
    output = subprocess.check_output(['git', 'show', revision_path])
    return json.loads(output), {'git_object': revision_path, 'bytes': len(output),
                                 'sha256': hashlib.sha256(output).hexdigest()}


def case_names(path, expected_count):
    cases = [test.attrib['name'] for test in ET.parse(path).iter('testcase')]
    require(len(cases) == len(set(cases)) == expected_count, f'Unexpected original roster in {path}')
    return cases


def elf_dependencies(executables):
    dependencies = set()
    for executable in executables:
        output = subprocess.check_output(['/usr/bin/ldd', str(executable)], text=True)
        require('not found' not in output, f'Unresolved dependency for {executable}')
        for line in output.splitlines():
            candidate = line.split('=>', 1)[-1].strip().split(' ', 1)[0]
            if candidate.startswith('/'):
                dependencies.add(Path(candidate).resolve(strict=True))
    return sorted(dependencies)


def verify_junit(path, expected):
    xml = ET.parse(path).getroot()
    for node in [xml, *xml.findall('.//testsuite')]:
        require(all(int(node.get(key, '0')) == 0 for key in CHECK_KEYS), f'Failed JUnit counters in {path}')
    cases = list(xml.iter('testcase'))
    require(Counter(c.attrib['name'] for c in cases) == Counter(expected), f'Unexpected JUnit roster in {path}')
    require(all(case.get('status') in (None, 'run') and
                not any(case.find(element) is not None for element in ('failure', 'error', 'skipped'))
                for case in cases), f'Incomplete or unsuccessful witness in {path}')
    return [{'name': case.attrib['name'], 'seconds': float(case.get('time', 0))} for case in cases]


def gtest_properties(directory, expected):
    records = {}
    for path in sorted(directory.glob('*.xml')):
        root = ET.parse(path).getroot()
        cases = list(root.iter('testcase'))
        require(len(cases) == 1, f'Expected one direct GTest witness in {path}')
        case, = cases
        name = case.attrib['classname'] + '.' + case.attrib['name']
        require(name not in records and case.get('status') == 'run' and case.get('result') == 'completed',
                f'Unexpected direct GTest record in {path}')
        require(all(int(node.get(k, '0')) == 0 for node in (root, *root.findall('testsuite'))
                    for k in CHECK_KEYS), f'GTest did not finish without failures or skips in {path}')
        records[name] = {'xml': pin(path), 'properties': {item.attrib['name']: item.attrib['value']
                                                        for item in case.findall('./properties/property')}}
    require(set(records) == set(expected), 'Original E2 producer XML is incomplete')
    return records


def collect_plan(model):
    criteria = {row['id']: row for row in model['acceptance_criteria']}
    dims = {row['id']: row for row in model['required_dimensions']}
    e2, e2_pin = git_object(E2_ACTIONS)
    names = [arg.removeprefix('--gtest_filter=') for action in e2.values()
             for arg in action['argv'] if arg.startswith('--gtest_filter=')]
    require(len(names) == 11 and set(names) == {row['name'] for row in criteria['E2']['evidence']},
            'E2 roster must contain all eleven original operating-model cases')
    authority = {'original_e2_actions': e2_pin}
    groups = {'original-e2': (names, 300, 4194304, False)}
    transport = {row['name'] for row in dims['polarised_transport_to_film']['evidence']}
    require(len(transport) == 10 and len(transport - set(names)) == 1,
            'Expected the original polarisation transport-to-film coverage')
    groups['e2-transport-supplement'] = (sorted(transport - set(names)), 300, 4194304, False)
    for label, size in (('original-53-cpu-cases', 53), ('original-connected-render', 1),
                        ('original-64x64-writers', 2)):
        historical = HISTORICAL / label
        names = case_names(historical / 'junit.xml', size)
        guardian = read(historical / 'report.json')
        require(guardian['minimum_host_available_kib'] == 2097152, f'Unexpected original reserve in {label}')
        authority[label] = [pin(historical / name) for name in ('junit.xml', 'report.json')]
        groups[label] = (names, guardian['maximum_seconds'], guardian['maximum_sampled_rss_kib'], False)

    p1 = {row['name'] for row in criteria['P1']['evidence']
          if row['name'].startswith('AnalyticValidationTest.')}
    p1.add('ShadowBoundary.KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p')
    p2 = {row['name'] for row in criteria['P2']['evidence']
          if row['name'].startswith(('BeamPropagationTest.', 'CpuJacobiOracle.'))}
    groups['p1-p2-host-oracles'] = (sorted(p1 | p2), 900, 4194304, False)

    p4, authority['original_p4_actions'] = git_object(P4_ACTIONS)
    original_p4 = next(iter(p4.values()))['argv'][-1]
    p4_cases = {row['name'] for group in (criteria['P4'], dims['page_thorne_and_volumetric_transfer'])
                for row in group['evidence'] if row['kind'] == 'gtest'}
    require(len([case for case in p4_cases if re.fullmatch(original_p4, case)]) == 3,
            'Expected all three original independent Page–Thorne controls')
    p4_device = {case for case in p4_cases if case.startswith('KernelParity.')}
    p4_cpu = {case for case in p4_cases - p4_device if not case.startswith('Vulkan')}
    groups['p4-host'] = (sorted(p4_cpu), 300, 4194304, False)
    groups['p4-kernel'] = (sorted(p4_device), 300, 4194304, True)

    output = {row['name'] for group in (criteria['P6'], dims['output_encoding'])
              for row in group['evidence'] if row['kind'] == 'gtest'}
    colour = {case for case in output if case.startswith('KernelParity.')}
    require(len(colour) == 2, 'Expected the original CIE and sRGB kernel cases')
    writers = set(groups['original-64x64-writers'][0])
    groups['p6-host'] = (sorted(output - colour - writers), 300, 4194304, False)
    groups['p6-kernel'] = (sorted(colour), 300, 4194304, True)
    order = ('original-e2', 'e2-transport-supplement', 'original-53-cpu-cases',
             'p1-p2-host-oracles', 'p4-host', 'p6-host', 'p6-kernel',
             'p4-kernel', 'original-connected-render', 'original-64x64-writers')
    return {key: groups[key] for key in order}, authority


def run():
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    require(not subprocess.check_output(['git', 'status', '--porcelain']), 'Science execution requires a clean tree')
    model = read(MODEL)
    rows = [row for area in ('acceptance_criteria', 'required_dimensions', 'capability_contracts')
            for row in model[area]]
    require(len(rows) == len({row['id'] for row in rows}) == 64, 'Operating-model row count or IDs drifted')
    groups, authority = collect_plan(model)
    work = WORK / 'runs'
    work.mkdir()
    RENDERS.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', work / 'run_owned.py')
    (work / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    discovery = json.loads(subprocess.check_output([*CTEST, '--show-only=json-v1']))
    records = {}
    for entry in discovery['tests']:
        records.setdefault(entry['name'], []).append(entry)
    all_cases = {case for roster, *_ in groups.values() for case in roster}
    registrations = {}
    for case in sorted(all_cases):
        options = records.get(case, [])
        require(len(options) == 1, f'Expected a unique registered case for {case}')
        entry, = options
        props = {item['name']: item['value'] for item in entry['properties']}
        require('Mandatory' in props.get('LABELS', []) and 'command' in entry, f'Case is not a built Mandatory test: {case}')
        registrations[case] = {'argv': entry['command'], 'labels': props['LABELS'],
                               'directory': props.get('WORKING_DIRECTORY')}
    write(work / 'registrations.json', registrations)
    executables = {Path(entry['argv'][0]) for entry in registrations.values()}
    runtime = BUILD / 'src/sirius/app'
    generated = BUILD / 'src/sirius/backend/retained'
    bound = sorted(executables | {runtime / 'sirius', MODEL, BUILD / 'CMakeCache.txt', BUILD / 'build.ninja',
                   BUILD / 'generated/sirius/alignment_receipt.json', work / 'registrations.json',
                   Path(__file__).resolve(), Path(CTEST[0]), Path('/usr/bin/cmake')})
    bound += sorted(path for place in (runtime / 'resources', BUILD / 'tests/backend/resources', generated)
                    for path in place.rglob('*') if path.is_file())
    bound += elf_dependencies(executables | {runtime / 'sirius'})
    bound = list(dict.fromkeys(bound))
    result = {'scope': 'Fresh original and operating-model science cohorts using development GCC products',
              'revision': revision, 'model': pin(MODEL), 'model_rows': len(rows),
              'authority': authority, 'registrations': pin(work / 'registrations.json'),
              'executable_inputs': [pin(path) for path in sorted(executables | {runtime / 'sirius'})],
              'build': {'configure_output': pin(WORK / 'configure.log'), 'build_output': pin(WORK / 'build.log'),
                        'starfield': pin(runtime / 'resources/assets/Starfield.png')},
              'cohorts': {}}
    result_path = WORK / 'result.json'
    actions = {}
    for label, (tests, seconds, memory, vulkan) in groups.items():
        selection = '^(' + '|'.join(re.escape(test) for test in tests) + ')$'
        exact = json.loads(subprocess.check_output([*CTEST, '-R', selection, '--show-only=json-v1']))
        require(Counter(entry['name'] for entry in exact['tests']) == Counter(tests),
                f'Exact CTest selection failed for {label}')
        directory = RENDERS / label
        directory.mkdir()
        argv = [*CTEST, '-V', '--no-tests=error', '--output-on-failure',
                '--output-junit', '{folder}/junit.xml', '-R', selection]
        if label == 'original-e2':
            argv = ['/usr/bin/cmake', '-E', 'env', 'GTEST_OUTPUT=xml:{folder}/gtest/', *argv]
        frozen = [*bound, *([ICD, PROVIDER, *elf_dependencies({PROVIDER})] if vulkan else [])]
        actions[label] = {'revision': revision, 'argv': argv, 'expected_exit': 0, 'kind': 'test',
                          'software_provider': vulkan, 'disable_software_shader_cache': True,
                          'runtime_provider': str(PROVIDER), 'runtime_icd': str(ICD),
                          'limits': {'seconds': seconds, 'rss_kib': memory, 'reserve_kib': 2097152},
                          'frozen_inputs': list(map(str, dict.fromkeys(frozen))),
                          'build_environment': {'TMPDIR': str(directory), 'MESA_SHADER_CACHE_DISABLE': 'true'},
                          'scope': 'Original or governed science cohort: ' + label}
    write(work / 'actions.json', actions)
    write(result_path, result)
    for label, (tests, *_rest) in groups.items():
        folder = work / label
        completed = subprocess.run([sys.executable, '-B', str(work / 'run_owned.py'), label])
        owner_path = folder / 'owner.json'
        owner = read(owner_path)
        cohort = {'owner': pin(owner_path), 'stdout': pin(folder / 'stdout.log'),
                  'stderr': pin(folder / 'stderr.log'),
                  'status': {key: owner[key] for key in ('passed', 'status', 'wall_seconds',
                             'peak_sampled_rss_kib', 'limits', 'stop_reason', 'child_exit', 'source_unchanged',
                             'whole_inputs_unchanged', 'cleanup_errors', 'remaining_owned_processes')}}
        junit = folder / 'junit.xml'
        if junit.is_file():
            cohort['junit'] = pin(junit)
        result['cohorts'][label] = cohort
        write(result_path, result)
        require(completed.returncode == 0 and owner['passed'] and owner['status'] == 'terminal' and
                owner['whole_inputs_unchanged'] and owner['source_unchanged'] and
                not owner['cleanup_errors'] and not owner['remaining_owned_processes'],
                f'Owned execution did not pass in {label}')
        cohort['cases'] = verify_junit(junit, tests)
        if label == 'original-e2':
            cohort['individual_witnesses'] = gtest_properties(folder / 'gtest', tests)
        if actions[label]['software_provider']:
            require(owner['provider_mappings'], f'Selected software Vulkan mapping missing in {label}')
            cohort['provider'] = pin(PROVIDER)
            cohort['icd'] = pin(ICD)
        cohort['accepted'] = True
        write(result_path, result)


if __name__ == '__main__':
    require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB available on /')
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    run()
