"""Exercise the original catalogue, governance, precision and disk cases on strict Linux."""

from collections import Counter
from pathlib import Path
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT / 'out/science-replay'))
from run_science import GUARD, ICD, MODEL, PROVIDER, elf_dependencies, pin, read, require, verify_junit, write

BUILD = ROOT / 'bin/linux-ci'
WORK = ROOT / 'out/strict-linux-science'
RENDERS = ROOT / 'renders/strict-linux-science'
CTEST = ['/usr/bin/ctest', '--test-dir', str(BUILD), '-C', 'Release']


def plan(model):
    criteria = {item['id']: item for item in model['acceptance_criteria']}
    e1 = {entry['name'] for entry in criteria['E1']['evidence']}
    e4 = {entry['name'] for entry in criteria['E4']['evidence']}
    e1_device = {name for name in e1 if name.startswith(('KernelParity.', 'VulkanRenderSession.'))}
    e4_device = {name for name in e4 if name.startswith('KernelParity.')}
    require(len(e1) == 28 and len(e4) == 7 and len(e1_device) == 2 and len(e4_device) == 1,
            'The original E1/E4 roster changed')
    p4_live = {entry['name'] for entry in criteria['P4']['evidence']
               if entry['name'].startswith('VulkanRenderSession.')}
    require(len(p4_live) == 2, 'The two original production disk witnesses are required')
    p1_host = {entry['name'] for entry in criteria['P1']['evidence']
               if entry['name'].startswith(('CpuTraceBoundary.', 'KerrTests.'))
               or entry['name'] == 'RK45IntegratorTests.AcceptedDenseErrorOnlyLimitsGrowth'}
    p1_device = {'KernelParity.PrecisionRungsConserveNearExtremalKerrWithoutImageComparison',
                 'KernelParity.PrecisionProbeArtifactsCarryOnlyTheirDeclaredFloat64Capability',
                 'VulkanBackend.KernelPrecisionDeclinesFmaBeforeDriverWorkUnlessBinary32WasEnabled',
                 'VulkanBackend.KernelPrecisionDeclinesUnsupportedFloat64AndMalformedInstructions'}
    require(p1_device <= {e['name'] for e in criteria['P1']['evidence']}, 'P1 witnesses must remain governed')
    sessions = {name for name in e1_device if name.startswith('VulkanRenderSession.')}
    turbulence = {name for name in p4_live if 'Turbulence' in name}
    supplemental = {
        'OperationalAttestation.FalseExternalEvidenceIsRejected',
        'OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice',
    }
    return {
        'source-and-release-governance': (e4 - e4_device | supplemental, False, 300),
        'e1-cpu-catalogue': (e1 - e1_device, False, 1800),
        'p1-host-boundaries': (p1_host, False, 1800),
        'p1-precision-admission': (p1_device - {'KernelParity.PrecisionRungsConserveNearExtremalKerrWithoutImageComparison'}, True, 900),
        'p1-three-rungs': ({'KernelParity.PrecisionRungsConserveNearExtremalKerrWithoutImageComparison'}, True, 3600),
        'e1-e4-kernel-oracles': ((e1_device - sessions) | e4_device, True, 1800),
        'e1-morris-thorne-parity': (sessions, True, 3600),
        'p4-live-turbulence': (turbulence, True, 3600),
        'p4-live-doppler': (p4_live - turbulence, True, 7200),
    }


def run():
    revision = subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip()
    require(not subprocess.check_output(['git','status','--porcelain']), 'Test from an exact clean tree')
    model = read(MODEL)
    require(sum(len(model[key]) for key in ('acceptance_criteria','required_dimensions','capability_contracts')) == 64,
            'Operating-model roster is not the original 64 rows')
    receipt = BUILD / 'generated/sirius/alignment_receipt.json'
    alignment = read(receipt)
    require(alignment['source_revision'] == revision and alignment['source_tree_clean'] is True and
            alignment['alignment_mode'] == 'qualification', 'Expected exact clean qualification candidate')
    work = WORK / 'runs'
    work.mkdir()
    RENDERS.mkdir()
    shaders = WORK / 'mesa-cache'
    shaders.mkdir()
    shutil.copyfile(GUARD/'run_owned.py', work/'run_owned.py')
    (work/'owned_guard.py').symlink_to(GUARD/'owned_guard.py')
    groups = plan(model)
    inventory = subprocess.check_output([*CTEST, '--show-only=json-v1'])
    (work/'inventory.json').write_bytes(inventory)
    import json
    rows = json.loads(inventory)['tests']
    require(len(rows) == len({row['name'] for row in rows}) and all('command' in row for row in rows),
            'Strict candidate registration is incomplete or duplicated')
    all_names = {name for cases, *_ in groups.values() for name in cases}
    selected = {}
    mandatory = []
    for row in rows:
        props = {entry['name']: entry['value'] for entry in row['properties']}
        if 'Mandatory' not in props.get('LABELS', []):
            continue
        mandatory.append(row['name'])
        if row['name'] in all_names:
            selected[row['name']] = {key: row[key] for key in ('name','command')}
            selected[row['name']]['working_directory'] = props.get('WORKING_DIRECTORY')
    require(set(selected) == all_names, f'Missing selected registrations: {sorted(all_names - set(selected))}')
    write(work/'selected.json', selected)
    runtime = BUILD / 'src/sirius/app'
    products = {runtime / name for name in ('sirius','sirius_app_tests','sirius_render_tests')}
    products |= {BUILD/'tests'/name/f'sirius_{name}_tests' for name in ('backend','core','oracle')}
    products.add(BUILD/'tests/sirius_base_tests')
    libraries = elf_dependencies(products)
    resources = runtime/'resources'
    require((resources/'assets/Starfield.png').samefile(ROOT/'assets/Starfield.png'), 'Use the existing starfield')
    static = sorted(products) + libraries + sorted(path for path in resources.rglob('*') if path.is_file())
    static += [MODEL, receipt, work/'inventory.json', work/'selected.json', Path(__file__).resolve(),
               BUILD/'build.ninja', BUILD/'CMakeCache.txt', Path(CTEST[0]), Path('/usr/bin/cmake')]
    static += sorted((BUILD/'src/sirius/backend/retained').glob('*.spv'))
    actions = {}
    for label, (roster, vulkan, seconds) in groups.items():
        path = RENDERS / label
        path.mkdir()
        regex = '^(' + '|'.join(re.escape(name) for name in sorted(roster)) + ')$'
        actual = json.loads(subprocess.check_output([*CTEST, '-R', regex, '--show-only=json-v1']))['tests']
        require(Counter(row['name'] for row in actual) == Counter(roster), f'Incorrect selection for {label}')
        actions[label] = {
            'revision': revision, 'kind': 'test', 'expected_exit': 0,
            'argv': [*CTEST, '-R', regex, '-V', '--output-on-failure', '--no-tests=error',
                     '--output-junit', '{folder}/junit.xml'],
            'limits': {'seconds': seconds, 'rss_kib': 8388608, 'reserve_kib': 2097152},
            'build_environment': {'TMPDIR': str(path), 'MESA_SHADER_CACHE_DIR': str(shaders)},
            'software_provider': vulkan, 'disable_software_shader_cache': False,
            'runtime_provider': str(PROVIDER), 'runtime_icd': str(ICD),
            'frozen_inputs': list(map(str, dict.fromkeys([
                *static, *([PROVIDER, ICD, *elf_dependencies({PROVIDER})] if vulkan else [])]))),
            'scope': 'Strict Linux qualification candidate, exact registered cases: ' + label,
        }
    write(work/'actions.json', actions)
    result = {'revision': revision, 'mode': 'qualification', 'scope': 'Selected strict Linux cases',
              'alignment_receipt': pin(receipt), 'model': pin(MODEL), 'model_rows': 64,
              'registered_cases': len(rows), 'registered_mandatory_cases': len(mandatory),
              'registered_case_inventory': pin(work/'inventory.json'), 'selected_registrations': selected,
              'products': [pin(path) for path in sorted(products)], 'actions': pin(work/'actions.json'),
              'runner': pin(Path(__file__).resolve()), 'cohorts': {}}
    report = WORK/'results.json'
    write(report, result)
    for label, (roster, vulkan, _) in groups.items():
        execution = work/label
        finished = subprocess.run([sys.executable, '-B', str(work/'run_owned.py'), label])
        owner_file = execution/'owner.json'
        owner = read(owner_file)
        cohort = {'owner': pin(owner_file), 'stdout': pin(execution/'stdout.log'),
                  'stderr': pin(execution/'stderr.log'), 'argv': owner['argv'],
                  'execution_environment': actions[label]['build_environment'],
                  'status': {field: owner[field] for field in ('status','passed','child_exit','wall_seconds',
                             'peak_sampled_rss_kib','limits','stop_reason','source_unchanged','whole_inputs_unchanged',
                             'remaining_owned_processes','cleanup_errors')}}
        result['cohorts'][label] = cohort
        require(owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
                and not owner['remaining_owned_processes'] and not owner['cleanup_errors'],
                f'Ownership or input identity failed in {label}')
        cohort['input_snapshot'] = pin(execution/'inputs-before.json')
        require(cohort['input_snapshot']['sha256'] == pin(execution/'inputs-after.json')['sha256'],
                f'Input snapshots differ in {label}')
        junit = execution/'junit.xml'
        if junit.exists():
            cohort['junit'] = pin(junit)
        if not finished.returncode:
            try:
                cohort['cases'] = verify_junit(junit, roster)
                if vulkan:
                    require(owner['provider_mappings'], f'No observed provider mapping in {label}')
                    cohort['provider'] = pin(PROVIDER)
                    cohort['icd'] = pin(ICD)
                if label == 'p1-three-rungs':
                    xml = execution/'gtest.xml'
                    properties = {item.attrib['name']: item.attrib['value']
                                  for item in ET.parse(xml).iter('property')}
                    cohort['precision_properties'] = properties
                    cohort['direct_test_xml'] = pin(xml)
                    require(properties.get('fp64_evidence') == 'native_invariants_executed',
                            'The three-rung witness did not execute native fp64')
                cohort['accepted'] = True
            except (RuntimeError, OSError, ET.ParseError) as error:
                cohort['validation_failure'] = str(error)
        if not cohort.get('accepted'):
            cohort['stdout_tail'] = (execution/'stdout.log').read_text().splitlines()[-20:]
            cohort['stderr_tail'] = (execution/'stderr.log').read_text().splitlines()[-10:]
        write(report, result)
    return all(record.get('accepted') for record in result['cohorts'].values())


if __name__ == '__main__':
    require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB available on /')
    sys.exit(0 if run() else 1)
