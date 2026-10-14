"""Seal three original finite P1/E4 CPU oracle controls without executing them."""
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
CASE = 'original-three-p1-e4-oracle-controls'
SOURCES = ['tests/oracle/symplectic_integrator_test.cpp', 'tests/oracle/polarisation_transport_test.cpp']
NAMES = sorted([
    'AnalyticValidationTest.NearExtremalKerrConservesEnergyAngularMomentumAndCarter',
    'AnalyticValidationTest.SymplecticChartExitAndNullProjectionFailureRemainDistinct',
    'OracleConnection.AnalyticConnectionAgreesWithMetricDerivatives',
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
    'prepare.py', 'run_oracle.py', 'run_owned.py', 'owned_guard.py', 'reference'}
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
                               'One frozen three-case original P1/E4 CPU oracle action')
wrapper = wrapper.replace('original-four-kerr-authority-controls', CASE).replace('"p4_complete": False', '"p1_complete": False, "e4_complete": False')
assert wrapper.count('len(NAMES) == 4') == wrapper.count('len(cases) == 4') == 1
wrapper = wrapper.replace('len(NAMES) == 4', 'len(NAMES) == 3').replace('len(cases) == 4', 'len(cases) == 3')
(WORK / 'run_oracle.py').write_text(wrapper)
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
criteria = [row for row in model['acceptance_criteria'] if row['id'] in ('P1', 'E4')]
assert set(NAMES) <= {row['name'] for criterion in criteria for row in criterion['evidence'] if row['kind'] == 'gtest'}
producer = dict(base['source_and_registration_join']['producers']['oracle'])
assert producer == {'path': 'bin/linux-gcc/tests/oracle/sirius_oracle_tests',
    'recovery_archive': 'attestations/alignment/4bcd987/gcc-revision-synchronisation/prior-versions/be9a6f6136937567e9537e4ff25ee71589c1f98c46aa04b25790a922832abba4',
    'bytes': 1678696, 'sha256': 'be9a6f6136937567e9537e4ff25ee71589c1f98c46aa04b25790a922832abba4'}
assert identity(ROOT / producer['path']) == identity(ROOT / producer['recovery_archive']) == {
    field: producer[field] for field in ('bytes', 'sha256')}
review_path = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14/review.json'
record = json.loads(review_path.read_text())['current_live_byte_bindings']['tested_artifacts']['records']['sirius_oracle_tests']
assert all(record[field] == producer[field] for field in ('bytes', 'sha256'))
accepted_path = ROOT / base['source_and_registration_join']['accepted_source_build_evidence']
accepted = json.loads(accepted_path.read_text())
assert accepted['passed'] and accepted['source_revision'] == PRODUCT and accepted['source_unchanged']
# Snapshot the original checked-in independent reference and its local include
# closure; no replacement oracle or new numerical threshold is introduced.
snapshot_names = set()
pending = [*SOURCES, 'tests/oracle/CMakeLists.txt', 'tests/operating_model.json',
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
compiles = []
source_for_name = {}
for source in SOURCES:
    selected = [row for row in commands if Path(row['file']).resolve() == ROOT / source]
    assert len(selected) == 1
    row = selected[0]
    assert row['directory'] == str(ROOT / 'bin/linux-gcc')
    assert row['output'] == 'tests/oracle/CMakeFiles/sirius_oracle_tests.dir/' + Path(source).name + '.o'
    for flag in ('/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-DSIRIUS_RELEASE_RESOURCE_LOCKED=1',
                 '-O3', '-DNDEBUG', '-std=c++2c', '-Werror', 'tests/support/fail_on_skip.h'):
        assert flag in row['command'], (source, flag)
    assert Path(source).name in (ROOT / 'tests/oracle/CMakeLists.txt').read_text()
    source_names = {'{}.{}'.format(*pair) for pair in re.findall(r'TEST\(\s*(\w+)\s*,\s*(\w+)\s*\)', (ROOT / source).read_text())}
    for name in set(NAMES) & source_names:
        assert name not in source_for_name
        source_for_name[name] = source
    compiles.append(row)
assert set(source_for_name) == set(NAMES)
registration = ROOT / 'bin/linux-gcc/tests/oracle/sirius_oracle_tests[1]_tests.cmake'
registrations = []
for name in NAMES:
    lines = registration.read_text().splitlines()
    adds = [line for line in lines if line.startswith('add_test([=[' + name + ']=]')]
    props = [line for line in lines if line.startswith('set_tests_properties([=[' + name + ']=]')]
    assert len(adds) == len(props) == 1 and str(ROOT / producer['path']) + ' ' in adds[0]
    assert '--gtest_filter=' + name + ']==]' in adds[0] and 'LABELS Correctness' in props[0] and 'TIMEOUT' not in props[0]
    registrations.append({'name': name, 'source': source_for_name[name], 'registration_file': key(registration),
                          'add_test': adds[0], 'properties': props[0], 'producer': producer['path']})
links = [line for line in (ROOT / 'bin/linux-gcc/build.ninja').read_text().splitlines()
         if line.startswith('build tests/oracle/sirius_oracle_tests ') and ': CXX_EXECUTABLE_LINKER__' in line]
assert len(links) == 1 and all(row['output'] in links[0] for row in compiles)
cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
for name, value in base['source_and_registration_join']['cache_required'].items():
    assert re.search(r'^' + re.escape(name) + r':[^=]+=' + re.escape(value) + '$', cache, re.M), name
ctest = re.search(r'^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$', cache, re.M)
assert ctest
scopes = {
    'AnalyticValidationTest.NearExtremalKerrConservesEnergyAngularMomentumAndCarter': {
        'criterion': 'P1', 'metric': {'mass': 1, 'spin': 0.998},
        'initial_position_t_r_theta_phi': [0, 20, 1.1, 0],
        'initial_covariant_momentum': {'p_t': -1, 'p_theta': 0.4, 'p_phi': 3,
            'p_r': 'positive sqrt(-non_radial / inverse[1][1]) from original KerrMetricD inverse'},
        'integration': 'Original adaptive Hamiltonian absolute and relative tolerances1e-13; up to250 complete unit-affine spans, stopping on r>=200. Original default step and maximum-step controls remain unchanged.',
        'assertions': 'Every unit span complete; final radius>=200; maximum absolute E drift<1e-10, relative Lz drift<1e-10 and Q drift/max(abs(initial_Q),1e-30)<1e-10.',
        'scope_boundary': 'One original outgoing near-extremal Boyer-Lindquist ray. Oracle-side invariant witness, not live Cartesian/Vulkan conservation, capture, 1080p shadow or a full P1 claim.'},
    'AnalyticValidationTest.SymplecticChartExitAndNullProjectionFailureRemainDistinct': {
        'criterion': 'P1',
        'fixtures': 'Original flat near-axis Yoshida4 step1e-3 with null projection disabled; Schwarzschild impossible-null-root, infinite radial momentum and invalid axis near horizon; explicit periodic-null-projection schedule.',
        'assertions': 'Chart exit remains distinct from constraint failure; integrate retains last represented state and zero accepted steps on chart exit; invalid chart state is not capture; identical explicit schedule state produces identical q/p; accepted_steps_before=1 produces original constraint failure.',
        'scope_boundary': 'Original finite outcome/refusal/history controls only; no added fixture, domain or numerical threshold.'},
    'OracleConnection.AnalyticConnectionAgreesWithMetricDerivatives': {
        'criterion': 'E4', 'metric': {'mass': 1, 'spin': 0.998},
        'events_t_r_theta_phi': [[0, 8, 1.2, 0.3], [0, 4, 'pi/2 -0.2', 1], [0, 15, 0.7, 0]],
        'reference': 'Original NumericChristoffel central differences of KerrMetricD.Evaluate with h=1e-6*max(1,abs(x[sigma]))+1e-8, contracted with the same sampled inverse metric.',
        'assertions': 'At each of3 original events, maximum absolute difference over64 Christoffel components from KerrChristoffel <1e-9.',
        'independence_limit': 'Independent metric differentiation and connection assembly use the same Boyer-Lindquist metric/inverse authority. No independent all-event metric oracle or E2 transport control is selected.'},
}
assert set(scopes) == set(NAMES)
historical = ROOT / 'attestations/retained-p1/eb539fc/local-diagnostics/identity.json'
old_xml = historical.parent / 'oracle-near-extremal.xml'
history = json.loads(historical.read_text())
assert history['oracle_tests_sha256'] == 'a0c82a309ae6f3005a60fd820c86bdb0e03924869cb5f330e54ed978fb171f4d'
assert history['oracle_tests_sha256'] != producer['sha256']
assert identity(old_xml)['sha256'] == history['artifacts'][old_xml.name]
for path in (historical, old_xml):
    target = WORK / 'reference/historical' / path.name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)
    pins[key(path)], pins[key(target)] = identity(path), identity(target)
extra = [ROOT / producer['path'], ROOT / producer['recovery_archive'], Path(ctest.group(1)),
         Path('/usr/bin/python3.12'), review_path, accepted_path, registration, OWNER / 'manifest.json',
         OWNER / 'result.json', OWNER_RAW / 'preparation.json', OWNER_RAW / 'run_owned.py',
         OWNER_RAW / 'owned_guard.py', wrapper_path,
         *[WORK / name for name in ('prepare.py', 'run_oracle.py', 'run_owned.py', 'owned_guard.py')]]
for path in extra:
    pins[key(path)] = identity(path)
action = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'scope': 'Three original finite P1/E4 CPU oracle controls: one near-extremal invariant ray, finite symplectic chart/null/history refusals and3 connection events; no E2 control, all-event proof, full P1/E4/science/Mandatory or release qualification.',
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
    'ledger_scope': {'criteria': ['P1', 'E4'], 'selected_original_names': NAMES, 'selected_count': 3,
        'matching_current_producer_controls_repeated': False, 'p1_complete_claimed': False, 'e4_complete_claimed': False},
    'reference': {'original_source_snapshots': snapshots,
        'independent_reference': 'Original NumericChristoffel central metric differences and connection assembly plus original Hamiltonian E/Lz/Carter invariant calculation. The symplectic control independently asserts represented outcome/state/history behavior. No new reference implementation or threshold.',
        'independence_limit': 'The connection reference shares KerrMetricD metric/inverse values with its analytic connection; conservation uses the original Hamiltonian-state conserved-quantity formulas. These finite oracle checks do not establish live CPU/device propagation or an all-event mathematical proof.',
        'scope_boundary': 'Only the original three named P1/E4 witnesses, fixtures, steps and tolerances. No E2 transport cases, frames or stronger gates.',
        'historical_execution': {'xml': key(old_xml), 'identity': key(historical),
            'source_revision': history['source_revision'], 'oracle_sha256': history['oracle_tests_sha256'],
            'historical_near_extremal_seconds': 0.002,
            'reuse_boundary': 'The one preserved historical near-extremal pass uses a different oracle ELF; it is retained as historical. Bounded XML suffix search found no matching current-producer execution of any selected name, and no preserved XML occurrence of the other two names.'}},
    'source_and_registration_join': {'test_sources': SOURCES, 'independent_reference_sources': ['tests/oracle/polarisation_transport_test.cpp', 'src/sirius/oracle/transport_types.h'],
        'production_headers': sorted(name for name in snapshots if name.startswith('src/')),
        'producers': {'oracle': producer},
        'accepted_source_build_evidence': key(accepted_path),
        'cache_required': base['source_and_registration_join']['cache_required'],
        'original_compile_commands': compiles, 'original_link_ninja_lines': links,
        'registrations': registrations, 'selected_core_cases': 0, 'selected_oracle_cases': 3,
        'original_test_TIMEOUT_added': False},
    'owner_reuse': {'archive': key(OWNER), 'manifest_sha256': identity(OWNER / 'manifest.json')['sha256'],
        **{name: {'source': key(OWNER_RAW / name), **identity(WORK / name)} for name in ('run_owned.py', 'owned_guard.py')}},
    'private_tmpdir': {'path': str(WORK / 'temporary'), 'created_only_at_execution': True, 'mode': '0700',
        'owner_uid': os.getuid(), 'preserve_natural_residuals': True}, 'pins': pins})
for path in WORK.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
print(json.dumps({'prepared_without_execution': True, 'names': NAMES, 'pins': len(pins),
    'producer': producer, 'preparation': identity(WORK / 'preparation.json'),
    'wrapper': identity(WORK / 'run_oracle.py'),
    'command': ['python3.12', '-B', key(WORK / 'run_oracle.py')]}, indent=2))
