"""Freeze the three approved original Page-Thorne controls; no build or execution."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PRODUCT = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
RUNNER = 'f5d649bf5067d1d68dc9215462f37830122f0540'
CASE = 'original-three-page-thorne-controls'
NAMES = sorted('AnalyticValidationTest.' + name for name in (
    'PageThorneFluxMatchesIndependentQuadrature',
    'TruncatedPageThorneDiskUsesDeclaredZeroTorqueEdge',
    'PageThorneFluxApproachesNewtonianCubicFalloff'))
OWNER = ROOT / 'attestations/scientific-core/4bcd987/p4-cpu-radiation-current-source'
OWNER_RAW = OWNER / 'files/out/p4-cpu-radiation-current-source'
SOURCE = 'tests/oracle/analytic_benchmarks_test.cpp'
HEADERS = ['src/sirius/core/' + name for name in (
    'disk/novikov_thorne_disk.h', 'disk/disk_model.h', 'constants.h', 'coordinates.h',
    'kerr_orbits.h', 'spectral/spectral_radiance.h')]
HEADERS += ['src/sirius/oracle/beam_integrator.h', 'src/sirius/oracle/kerr_boyer_lindquist.h']


def identity(path):
    raw = Path(path).read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, timeout=30)


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + '\n')


assert Path.cwd() == ROOT and set(p.name for p in WORK.iterdir()) == {'prepare.py'}
assert git('rev-parse', 'HEAD').decode().strip() == RUNNER and not git('status', '--porcelain').strip()
changed = git('diff', '--name-only', PRODUCT, RUNNER).decode().splitlines()
assert set(changed) == {'.github/workflows/ci.yml', 'scripts/diagnose-macos-retained.py',
    'scripts/governance/build_policy.py', 'scripts/macos_metal_probe.py',
    'scripts/macos_moltenvk_integration.py', 'scripts/macos_vulkan_identity.py'}
assert identity(OWNER / 'manifest.json')['sha256'] == '982dbd299ab2ae9ef6ff5f9fd88d7c93efd0aef26e90c009b4392026479ea872'
records = {row['path']: {field: row[field] for field in ('bytes', 'sha256')}
           for row in json.loads((OWNER / 'manifest.json').read_text())['files']}
base = json.loads((OWNER_RAW / 'preparation.json').read_text())
assert identity(OWNER_RAW / 'preparation.json') == records[key(OWNER_RAW / 'preparation.json').removeprefix(key(OWNER) + '/')]
pins = {name: pin for name, pin in base['pins'].items()
        if name.startswith(('bin/', '/usr/', 'attestations/alignment/', 'attestations/fresh-session-review/'))}
for name, pin in pins.items():
    assert identity(ROOT / name) == pin, name
for name in ('run_owned.py', 'owned_guard.py'):
    source = OWNER_RAW / name
    assert identity(source) == records[str(source.relative_to(OWNER))]
    shutil.copy2(source, WORK / name)
wrapper_source = (OWNER_RAW / 'run_radiation.py').read_text()
tree = ast.parse(wrapper_source)
names_node = next(node for node in tree.body if isinstance(node, ast.Assign) and
                  any(isinstance(target, ast.Name) and target.id == 'NAMES' for target in node.targets))
lines = wrapper_source.splitlines(keepends=True)
lines[names_node.lineno - 1:names_node.end_lineno] = ['NAMES = ' + repr(set(NAMES)) + '\n']
wrapper = ''.join(lines).replace('One frozen nine-case original CPU radiation action', 'One frozen three-case original Page-Thorne action')
wrapper = wrapper.replace('original-nine-cpu-radiation-controls', CASE)
assert wrapper.count('len(NAMES) == 9') == wrapper.count('len(cases) == 9') == 1
wrapper = wrapper.replace('len(NAMES) == 9', 'len(NAMES) == 3').replace('len(cases) == 9', 'len(cases) == 3')
ast.parse(wrapper, filename='run_page_thorne.py')
(WORK / 'run_page_thorne.py').write_text(wrapper)
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
p4 = next(row for row in model['acceptance_criteria'] if row['id'] == 'P4')
assert set(NAMES) <= {row['name'] for row in p4['evidence'] if row['kind'] == 'gtest'}
producer = ROOT / 'bin/linux-gcc/tests/oracle/sirius_oracle_tests'
producer_pin = {'bytes': 1678696, 'sha256': 'be9a6f6136937567e9537e4ff25ee71589c1f98c46aa04b25790a922832abba4'}
archive = ROOT / 'attestations/alignment/4bcd987/gcc-revision-synchronisation/prior-versions' / producer_pin['sha256']
assert identity(producer) == identity(archive) == producer_pin
review_path = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14/review.json'
review = json.loads(review_path.read_text())
assert {name: review['current_live_byte_bindings']['tested_artifacts']['records']['sirius_oracle_tests'][name]
        for name in ('bytes', 'sha256')} == producer_pin
snapshots = {}
for name in [SOURCE, *HEADERS, 'tests/oracle/CMakeLists.txt', 'tests/operating_model.json',
             'tests/labels/CTestLabels.cmake', 'tests/support/fail_on_skip.h']:
    path = ROOT / name
    assert git('show', PRODUCT + ':' + name) == path.read_bytes(), name
    snapshots[name] = {**identity(path), 'git_blob': git('rev-parse', PRODUCT + ':' + name).decode().strip()}
    target = WORK / 'reference' / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)
    pins[name], pins[key(target)] = identity(path), identity(target)
commands = json.loads((ROOT / 'bin/linux-gcc/compile_commands.json').read_text())
compiles = [row for row in commands if Path(row['file']).resolve() == ROOT / SOURCE]
assert len(compiles) == 1
compile_row = compiles[0]
assert compile_row['directory'] == str(ROOT / 'bin/linux-gcc')
assert 'sirius_oracle_tests.dir/' in compile_row['output'] and compile_row['output'].endswith('analytic_benchmarks_test.cpp.o')
for flag in ('/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-O3', '-DNDEBUG', '-std=c++2c', '-Werror', 'tests/support/fail_on_skip.h'):
    assert flag in compile_row['command'], flag
source_names = {'{}.{}'.format(*pair) for pair in re.findall(r'TEST\(\s*(\w+)\s*,\s*(\w+)\s*\)', (ROOT / SOURCE).read_text())}
assert set(NAMES) <= source_names
registration = ROOT / 'bin/linux-gcc/tests/oracle/sirius_oracle_tests[1]_tests.cmake'
lines = registration.read_text().splitlines()
registrations = []
for name in NAMES:
    adds = [line for line in lines if line.startswith('add_test([=[' + name + ']=]')]
    props = [line for line in lines if line.startswith('set_tests_properties([=[' + name + ']=]')]
    assert len(adds) == len(props) == 1 and str(producer) + ' ' in adds[0]
    assert '--gtest_filter=' + name + ']==]' in adds[0] and 'LABELS Correctness' in props[0] and 'TIMEOUT' not in props[0]
    registrations.append({'name': name, 'source': SOURCE, 'registration_file': key(registration),
                          'add_test': adds[0], 'properties': props[0], 'producer': key(producer)})
links = [line for line in (ROOT / 'bin/linux-gcc/build.ninja').read_text().splitlines()
         if line.startswith('build tests/oracle/sirius_oracle_tests ') and ': CXX_EXECUTABLE_LINKER__' in line]
assert len(links) == 1 and compile_row['output'] in links[0]
cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
for name, value in base['source_and_registration_join']['cache_required'].items():
    assert re.search(r'^' + re.escape(name) + r':[^=]+=' + re.escape(value) + '$', cache, re.M), name
ctest = re.search(r'^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$', cache, re.M)
assert ctest
scopes = {
    'AnalyticValidationTest.PageThorneFluxMatchesIndependentQuadrature': {'spins': [-0.7, 0.0, 0.9], 'radius_over_isco': [1.05, 1.2, 2.0, 4.0],
        'reference_radius_over_isco': 1.5, 'sample_pairs': 12,
        'assertions': 'Positive reference temperature; original EXPECT_NEAR(production/oracle,1,2e-5) for normalized temperatures.',
        'independence_limit': 'Reference inner edge uses production ComputeIsco(spin); independence concerns flux integration conditional on that edge, not an independent ISCO oracle.'},
    'AnalyticValidationTest.TruncatedPageThorneDiskUsesDeclaredZeroTorqueEdge': {
        'spin': 0.0, 'inner_radius': 10.0, 'outer_radius': 100.0, 'reference_radius': 15.0,
        'radius_over_inner': [1.05, 1.2, 2.0, 4.0],
        'assertions': 'Original exact InnerRadius=10 and zero Flux at10 and10*1.001; positive reference temperature; normalized-temperature ratio within2e-5 against independently integrated flux with literal inner edge10.'},
    'AnalyticValidationTest.PageThorneFluxApproachesNewtonianCubicFalloff': {
        'outer_radius': 100000.0, 'flux_radii': [10000.0, 20000.0],
        'assertions': 'Both fluxes positive; original EXPECT_NEAR(inner_flux/outer_flux,8,0.2); finite two-radius tail check only.'}}
assert set(scopes) == set(NAMES)
historical = ROOT / 'attestations/scientific-core/5e70607/final-build-identity-5e70607.json'
historical_identity = json.loads(historical.read_text())
assert historical_identity['artifacts'][key(producer)] != producer_pin['sha256']
old_xml = historical.parent / 'final-integrated-cpu.xml'
old_passes = [row.attrib for row in ET.parse(old_xml).getroot().iter('testcase') if row.get('name') in NAMES and
              not any(row.findall(tag) for tag in ('failure', 'error', 'skipped'))]
extra = [producer, archive, registration, Path(ctest.group(1)), Path('/usr/bin/python3.12'),
         review_path, OWNER / 'manifest.json', OWNER / 'disposition.json', OWNER_RAW / 'preparation.json',
         OWNER_RAW / 'run_owned.py', OWNER_RAW / 'owned_guard.py', OWNER_RAW / 'run_radiation.py', historical, old_xml,
         *[WORK / name for name in ('prepare.py', 'run_page_thorne.py', 'run_owned.py', 'owned_guard.py')]]
for path in extra:
    pins[key(path)] = identity(path)
action = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'scope': 'Three original finite CPU Page-Thorne oracle controls only; no device/frame, independent ISCO claim, full P4/science/Mandatory or release qualification.',
    'argv': [ctest.group(1), '--test-dir', 'bin/linux-gcc', '-V', '--output-on-failure', '--no-tests=error',
             '--output-junit', '{folder}/ctest.xml', '--output-log', '{folder}/ctest.log', '-R',
             '^(' + '|'.join(re.escape(name) for name in NAMES) + ')$'],
    'frozen_inputs': list(pins) + [key(WORK / 'actions.json'), key(WORK / 'preparation.json')]}
write('actions.json', {CASE: action})
pins[key(WORK / 'actions.json')] = identity(WORK / 'actions.json')
write('preparation.json', {
    'prepared_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'product_source_revision': PRODUCT,
    'runner_source_revision': RUNNER, 'diagnostic_only_difference': changed, 'executed': False,
    'expected_names': NAMES, 'scope': action['scope'], 'finite_acceptance_scope': scopes,
    'ledger_scope': {'criterion': 'P4', 'selected_original_names': NAMES, 'selected_count': 3, 'other_completed_controls_repeated': False},
    'reference': {'original_source_snapshots': snapshots,
        'independent_flux_reference': 'Original checked-in Page & Thorne1974 Eqs11b/15n flux shape q(r)/r,32768-panel midpoint quadrature with independent circular E,L,Omega and five-point finite differences at h=1e-4*r; production uses16-point Gauss-Legendre and analytic derivatives.',
        'scope_boundary': 'Finite normalized-temperature comparisons; shared production ISCO in the default-edge helper, literal inner edge10 in the truncated helper. No new reference implementation, threshold, all-spin/radius proof or absolute physical luminosity claim.',
        'historical_execution': {'xml': key(old_xml), 'identity': key(historical), 'revision': historical_identity['revision'],
            'oracle_producer_sha256': historical_identity['artifacts'][key(producer)], 'selected_old_passes': old_passes,
            'reuse_boundary': 'Two selected names pass historically on a different oracle ELF; the truncated-edge name is absent from this historical XML. No matching-current selected execution found in the bounded reviewed bundles; no historical result relabelled.'}},
    'source_and_registration_join': {'test_sources': [SOURCE], 'production_headers': HEADERS,
        'oracle_producer': producer_pin, 'oracle_recovery_archive': key(archive),
        'accepted_source_build_evidence': base['source_and_registration_join']['accepted_source_build_evidence'],
        'cache_required': base['source_and_registration_join']['cache_required'],
        'original_compile_commands': compiles, 'original_link_ninja_lines': links, 'registrations': registrations,
        'selected_oracle_cases': 3, 'original_test_TIMEOUT_added': False},
    'owner_reuse': {'archive': key(OWNER), 'manifest_sha256': identity(OWNER / 'manifest.json')['sha256'],
        **{name: {'source': key(OWNER_RAW / name), **identity(WORK / name)} for name in ('run_owned.py', 'owned_guard.py')}},
    'private_tmpdir': {'path': str(WORK / 'temporary'), 'created_only_at_execution': True, 'mode': '0700',
                      'owner_uid': os.getuid(), 'preserve_natural_residuals': True}, 'pins': pins})
for path in WORK.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
print(json.dumps({'prepared_without_execution': True, 'names': NAMES, 'pins': len(pins),
    'producer': producer_pin, 'preparation': identity(WORK / 'preparation.json'),
    'command': ['python3.12', '-B', key(WORK / 'run_page_thorne.py')]}, indent=2))
