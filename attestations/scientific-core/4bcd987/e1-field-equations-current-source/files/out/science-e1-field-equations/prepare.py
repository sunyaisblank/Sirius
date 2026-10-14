"""Seal two original finite E1 field-equation controls without executing them."""
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

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PRODUCT = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
RUNNER = 'f5d649bf5067d1d68dc9215462f37830122f0540'
CASE = 'original-two-field-equation-controls'
SOURCE = 'tests/core/kerr_schild_field_equation_test.cpp'
ORACLE = 'tests/core/metric_curvature_oracle.h'
NAMES = sorted([
    'KerrSchildFieldEquations.LiveCartesianFamilySatisfiesEinsteinMaxwell',
    'KerrSchildFieldEquations.KerrNewmanPotentialIsSourceFreeOutsideTheRing',
])
OWNER = ROOT / 'attestations/scientific-core/4bcd987/p4-kerr-authority-current-source'
OWNER_RAW = OWNER / 'files/out/science-p4-kerr-authority'


def identity(path):
    data = Path(path).read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, timeout=30)


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + '\n')


assert not sys.flags.optimize and not os.environ.get('PYTHONOPTIMIZE')
assert Path.cwd() == ROOT and {p.name for p in WORK.iterdir()} <= {
    'prepare.py', 'run_field_equations.py', 'run_owned.py', 'owned_guard.py', 'reference'}
assert git('rev-parse', 'HEAD').decode().strip() == RUNNER and not git('status', '--porcelain').strip()
changed = git('diff', '--name-only', PRODUCT, RUNNER).decode().splitlines()
assert set(changed) == {'.github/workflows/ci.yml', 'scripts/diagnose-macos-retained.py',
    'scripts/governance/build_policy.py', 'scripts/macos_metal_probe.py',
    'scripts/macos_moltenvk_integration.py', 'scripts/macos_vulkan_identity.py'}
assert identity(OWNER / 'manifest.json')['sha256'] == '137f1410f1087f213b2a9af58cd8c7fa5750ff5453e597be1e619a970d49258b'
records = json.loads((OWNER / 'manifest.json').read_text())['files']
assert identity(OWNER_RAW / 'preparation.json') == records[str((OWNER_RAW / 'preparation.json').relative_to(OWNER))]
base = json.loads((OWNER_RAW / 'preparation.json').read_text())
# Preserve existing CTest discovery and strict build/producer authorities, without
# importing unrelated frame, shader or asset inputs into this CPU-only action.
pins = {name: pin for name, pin in base['pins'].items()
        if name.startswith(('/usr/', 'attestations/fresh-session-review/'))
        or name == base['source_and_registration_join']['accepted_source_build_evidence']
        or name in ('bin/linux-gcc/CMakeCache.txt', 'bin/linux-gcc/build.ninja', 'bin/linux-gcc/compile_commands.json')
        or name.startswith('bin/linux-gcc/') and name.endswith('.cmake')}
for name, pin in pins.items():
    assert identity(ROOT / name) == pin, name
for name in ('run_owned.py', 'owned_guard.py'):
    source = OWNER_RAW / name
    assert identity(source) == records[str(source.relative_to(OWNER))]
    shutil.copy2(source, WORK / name)
wrapper_path = OWNER_RAW / 'run_kerr_authority.py'
assert identity(wrapper_path) == records[str(wrapper_path.relative_to(OWNER))]
wrapper = wrapper_path.read_text()
node = next(n for n in ast.parse(wrapper).body if isinstance(n, ast.Assign) and
            any(isinstance(t, ast.Name) and t.id == 'NAMES' for t in n.targets))
lines = wrapper.splitlines(keepends=True)
lines[node.lineno - 1:node.end_lineno] = ['NAMES = ' + repr(set(NAMES)) + '\n']
wrapper = ''.join(lines).replace('One frozen four-case original P4 Kerr authority action',
                               'One frozen two-case original E1 field-equation action')
wrapper = wrapper.replace('original-four-kerr-authority-controls', CASE).replace('"p4_complete"', '"e1_complete"')
assert wrapper.count('len(NAMES) == 4') == wrapper.count('len(cases) == 4') == 1
wrapper = wrapper.replace('len(NAMES) == 4', 'len(NAMES) == 2').replace('len(cases) == 4', 'len(cases) == 2')
(WORK / 'run_field_equations.py').write_text(wrapper)
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
e1 = next(row for row in model['acceptance_criteria'] if row['id'] == 'E1')
assert set(NAMES) <= {row['name'] for row in e1['evidence'] if row['kind'] == 'gtest'}
producer = dict(base['source_and_registration_join']['producers']['core'])
assert producer == {'path': 'bin/linux-gcc/tests/core/sirius_core_tests',
    'recovery_archive': 'attestations/alignment/4bcd987/gcc-revision-synchronisation/prior-versions/c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58',
    'bytes': 5199120, 'sha256': 'c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58'}
assert identity(ROOT / producer['path']) == identity(ROOT / producer['recovery_archive']) == {
    field: producer[field] for field in ('bytes', 'sha256')}
review_path = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14/review.json'
record = json.loads(review_path.read_text())['current_live_byte_bindings']['tested_artifacts']['records']['sirius_core_tests']
assert all(record[field] == producer[field] for field in ('bytes', 'sha256'))
accepted_path = ROOT / base['source_and_registration_join']['accepted_source_build_evidence']
accepted = json.loads(accepted_path.read_text())
assert accepted['passed'] and accepted['source_revision'] == PRODUCT and accepted['source_unchanged']
# Snapshot the original checked-in independent reference and its local include
# closure; no replacement oracle or new numerical threshold is introduced.
snapshot_names = set()
pending = [SOURCE, ORACLE, 'tests/core/CMakeLists.txt', 'tests/operating_model.json',
           'tests/labels/CTestLabels.cmake', 'tests/support/fail_on_skip.h']
while pending:
    name = pending.pop()
    if name in snapshot_names:
        continue
    snapshot_names.add(name)
    path = ROOT / name
    for include in re.findall(r'^#include\s+"([^"]+)"', path.read_text(), re.M):
        candidate = ROOT / 'src' / include if include.startswith('sirius/') else path.parent / include
        assert candidate.is_file(), (name, include)
        pending.append(key(candidate))
snapshots = {}
for name in sorted(snapshot_names):
    path = ROOT / name
    assert git('show', PRODUCT + ':' + name) == path.read_bytes(), name
    snapshots[name] = {**identity(path), 'git_blob': git('rev-parse', PRODUCT + ':' + name).decode().strip()}
    target = WORK / 'reference' / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)
    pins[name], pins[key(target)] = identity(path), identity(target)
commands = json.loads((ROOT / 'bin/linux-gcc/compile_commands.json').read_text())
selected = [row for row in commands if Path(row['file']).resolve() == ROOT / SOURCE]
assert len(selected) == 1
compile_row = selected[0]
assert compile_row['directory'] == str(ROOT / 'bin/linux-gcc')
assert compile_row['output'] == 'tests/core/CMakeFiles/sirius_core_tests.dir/kerr_schild_field_equation_test.cpp.o'
for flag in ('/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-DSIRIUS_RELEASE_RESOURCE_LOCKED=1',
             '-O3', '-DNDEBUG', '-std=c++2c', '-Werror', 'tests/support/fail_on_skip.h'):
    assert flag in compile_row['command'], flag
assert 'kerr_schild_field_equation_test.cpp' in (ROOT / 'tests/core/CMakeLists.txt').read_text()
source_names = {'{}.{}'.format(*pair) for pair in re.findall(r'TEST\(\s*(\w+)\s*,\s*(\w+)\s*\)', (ROOT / SOURCE).read_text())}
assert set(NAMES) == source_names
registration = ROOT / 'bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake'
registrations = []
for name in NAMES:
    lines = registration.read_text().splitlines()
    adds = [line for line in lines if line.startswith('add_test([=[' + name + ']=]')]
    props = [line for line in lines if line.startswith('set_tests_properties([=[' + name + ']=]')]
    assert len(adds) == len(props) == 1 and str(ROOT / producer['path']) + ' ' in adds[0]
    assert '--gtest_filter=' + name + ']==]' in adds[0] and 'LABELS Correctness' in props[0] and 'TIMEOUT' not in props[0]
    registrations.append({'name': name, 'source': SOURCE, 'registration_file': key(registration),
                          'add_test': adds[0], 'properties': props[0], 'producer': producer['path']})
links = [line for line in (ROOT / 'bin/linux-gcc/build.ninja').read_text().splitlines()
         if line.startswith('build tests/core/sirius_core_tests ') and ': CXX_EXECUTABLE_LINKER__' in line]
assert len(links) == 1 and compile_row['output'] in links[0]
cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
for name, value in base['source_and_registration_join']['cache_required'].items():
    assert re.search(r'^' + re.escape(name) + r':[^=]+=' + re.escape(value) + '$', cache, re.M), name
ctest = re.search(r'^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$', cache, re.M)
assert ctest
fixture_rows = re.findall(r'\{"([^"]+)",\s*(KerrSchildParams::\w+\([^)]*\)|KerrSchildParams\{[^}]+\}),\s*([\d.]+),\s*([\d.]+),\s*(-?[\d.]+)(?:,\s*(true|false),\s*(true|false))?\s*\}', (ROOT / SOURCE).read_text())
assert len(fixture_rows) == 14
fixtures = []
for label, constructor, radius, theta, phi, extremal, horizonless in fixture_rows:
    charged = 'ReissnerNordstrom' in constructor or 'KerrNewman' in constructor
    fixtures.append({'name': label, 'original_parameter_constructor': constructor,
        'radius_theta_phi': [float(radius), float(theta), float(phi)], 'charged': charged,
        'require_extremal': extremal == 'true', 'require_horizonless': horizonless == 'true'})
assert sum(row['charged'] for row in fixtures) == 8
scopes = {
    'KerrSchildFieldEquations.LiveCartesianFamilySatisfiesEinsteinMaxwell': {
        'fixture_count': 14, 'fixtures': fixtures,
        'curvature_step': '2e-4 * radius', 'potential_field_step': '5e-5 * radius',
        'tolerance': 'max(2e-8 / radius^2, 3e-4 * equation_scale)',
        'assertions': 'Original scalar R=4*Lambda and all16 Einstein-Maxwell components; original charged/cosmological source magnitude >10*tolerance guards. Original marked horizonless/extremal checks, extremality tolerance2e-15.',
        'domain_note': 'Stored binary64 a=.8,Q=.6 fixture is horizonless; exactly extremal fixture is M5,a4,Q3. No all-event or continuous-domain proof.'},
    'KerrSchildFieldEquations.KerrNewmanPotentialIsSourceFreeOutsideTheRing': {
        'fixture_count': 8, 'fixtures': [row['name'] for row in fixtures if row['charged']],
        'charged_subset': 'Original loop skips Q==0: only the eight charged fixtures, not all fourteen geometries.',
        'determinant_assertion': 'det(g) near -1 with2e-12', 'derivative_step': '2e-4 * radius',
        'potential_field_step': '4e-5 * radius', 'tolerance': 'max(2e-8, 2e-4 * abs(Q) / radius^3)',
        'assertions': 'Four original source-free divergence components of sqrt(-g)*F^{mu nu}; spatial coordinate derivatives use the original fourth-order centred stencil and stationarity.'},
}
extra = [ROOT / producer['path'], ROOT / producer['recovery_archive'], Path(ctest.group(1)),
         Path('/usr/bin/python3.12'), review_path, accepted_path, registration, OWNER / 'manifest.json',
         OWNER / 'result.json', OWNER_RAW / 'preparation.json', OWNER_RAW / 'run_owned.py',
         OWNER_RAW / 'owned_guard.py', wrapper_path,
         *[WORK / name for name in ('prepare.py', 'run_field_equations.py', 'run_owned.py', 'owned_guard.py')]]
for path in extra:
    pins[key(path)] = identity(path)
action = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'scope': 'Two original finite E1 CPU Cartesian Kerr-Schild field-equation controls:14 Einstein-Maxwell fixtures and8 charged Maxwell fixtures only; no all-event proof, full E1/science/Mandatory or release qualification.',
    'argv': [ctest.group(1), '--test-dir', 'bin/linux-gcc', '-V', '--output-on-failure', '--no-tests=error',
             '--output-junit', '{folder}/ctest.xml', '--output-log', '{folder}/ctest.log', '-R',
             '^(' + '|'.join(re.escape(name) for name in NAMES) + ')$'],
    'frozen_inputs': list(pins) + [key(WORK / 'actions.json'), key(WORK / 'preparation.json')]}
write('actions.json', {CASE: action})
pins[key(WORK / 'actions.json')] = identity(WORK / 'actions.json')
write('preparation.json', {'prepared_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'product_source_revision': PRODUCT, 'runner_source_revision': RUNNER,
    'diagnostic_only_difference': changed, 'executed': False, 'expected_names': NAMES,
    'scope': action['scope'], 'finite_acceptance_scope': scopes,
    'ledger_scope': {'criterion': 'E1', 'selected_original_names': NAMES, 'selected_count': 2,
        'other_completed_controls_repeated': False, 'e1_complete_claimed': False},
    'reference': {'original_source_snapshots': snapshots,
        'independent_reference': 'Original metric_curvature_oracle.h constructs the Levi-Civita connection from live g/dg and differentiates it with a fourth-order centred stencil. Original test separately differentiates A_mu=-Q*r*l_mu/Sigma, assembles F and8*pi*T, and checks its divergence; it never reuses production H or closed-form Ricci/stress.',
        'independence_limit': 'Live analytic metric/first derivatives enter the connection reference. Both sides use live sampled metric and generic TensorOps inversion to raise/contract; this is finite independent differential consistency, not an independent implementation of all geometry.',
        'scope_boundary': 'Exactly the original14 events and8 charged-event subset, stencils and tolerances. No new reference implementation, refinement, fixture, threshold or stronger gate.',
        'execution_gap': 'Bounded prior review found no matching preserved current-core execution of these two original names. This preparation executes nothing and does not relabel historical or other witnesses.'},
    'source_and_registration_join': {'test_sources': [SOURCE], 'independent_reference_sources': [ORACLE],
        'production_headers': sorted(name for name in snapshots if name.startswith('src/')),
        'producers': {'core': producer},
        'accepted_source_build_evidence': key(accepted_path),
        'cache_required': base['source_and_registration_join']['cache_required'],
        'original_compile_commands': [compile_row], 'original_link_ninja_lines': links,
        'registrations': registrations, 'selected_core_cases': 2, 'selected_oracle_cases': 0,
        'original_test_TIMEOUT_added': False},
    'owner_reuse': {'archive': key(OWNER), 'manifest_sha256': identity(OWNER / 'manifest.json')['sha256'],
        **{name: {'source': key(OWNER_RAW / name), **identity(WORK / name)} for name in ('run_owned.py', 'owned_guard.py')}},
    'private_tmpdir': {'path': str(WORK / 'temporary'), 'created_only_at_execution': True, 'mode': '0700',
        'owner_uid': os.getuid(), 'preserve_natural_residuals': True}, 'pins': pins})
for path in WORK.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
print(json.dumps({'prepared_without_execution': True, 'names': NAMES, 'pins': len(pins),
    'producer': producer, 'preparation': identity(WORK / 'preparation.json'),
    'wrapper': identity(WORK / 'run_field_equations.py'),
    'command': ['python3.12', '-B', key(WORK / 'run_field_equations.py')]}, indent=2))
