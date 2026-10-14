"""Seal four original P4 CPU orbit/ZAMO controls without executing them."""
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
CASE = 'original-four-kerr-authority-controls'
SOURCES = {
    'tests/core/kerr_orbit_authority_test.cpp': 'core',
    'tests/oracle/kerr_metric_test.cpp': 'oracle',
}
NAMES = sorted([
    'KerrOrbitAuthority.EveryCpuConsumerSharesSignedIscoAndCircularEmitterLaw',
    'KerrOrbitAuthority.UnrepresentedInputsDeclineInsteadOfChangingTheOrbit',
    'KerrZamoTransfer.ArbitraryLatitudeMatchesIndependentBoyerLindquistContraction',
    'KerrZamoTransfer.KerrSchildSlicingNormalIsNotTheOffEquatorialZamo',
])
HEADERS = ['src/sirius/' + name for name in (
    'core/kerr_orbits.h', 'core/relativistic_transfer.h', 'core/metrics/kerr_schild_family.h',
    'core/disk/novikov_thorne_disk.h', 'core/disk/disk_model.h', 'core/coordinates.h',
    'core/constants.h', 'core/tensor.h', 'core/dual_number.h', 'core/first_order_number.h',
    'core/second_order_number.h', 'core/twofold.h', 'core/metrics/metric.h',
    'core/metrics/registry.h', 'core/spectral/spectral_radiance.h',
    'oracle/kerr_boyer_lindquist.h', 'oracle/metric_interface.h', 'oracle/transport_types.h',
    'base/contracts.h')]
OWNER = ROOT / 'attestations/scientific-core/4bcd987/p4-cpu-radiation-current-source'
OWNER_RAW = OWNER / 'files/out/p4-cpu-radiation-current-source'


def identity(path):
    data = Path(path).read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, timeout=30)


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + '\n')


assert Path.cwd() == ROOT and {p.name for p in WORK.iterdir()} == {'prepare.py'}
assert git('rev-parse', 'HEAD').decode().strip() == RUNNER and not git('status', '--porcelain').strip()
changed = git('diff', '--name-only', PRODUCT, RUNNER).decode().splitlines()
assert set(changed) == {'.github/workflows/ci.yml', 'scripts/diagnose-macos-retained.py',
    'scripts/governance/build_policy.py', 'scripts/macos_metal_probe.py',
    'scripts/macos_moltenvk_integration.py', 'scripts/macos_vulkan_identity.py'}
assert identity(OWNER / 'manifest.json')['sha256'] == '982dbd299ab2ae9ef6ff5f9fd88d7c93efd0aef26e90c009b4392026479ea872'
records = {row['path']: {field: row[field] for field in ('bytes', 'sha256')}
           for row in json.loads((OWNER / 'manifest.json').read_text())['files']}
assert identity(OWNER_RAW / 'preparation.json') == records[str((OWNER_RAW / 'preparation.json').relative_to(OWNER))]
base = json.loads((OWNER_RAW / 'preparation.json').read_text())
pins = {name: pin for name, pin in base['pins'].items()
        if name.startswith(('bin/', '/usr/', 'attestations/alignment/', 'attestations/fresh-session-review/'))}
for name, pin in pins.items():
    assert identity(ROOT / name) == pin, name
for name in ('run_owned.py', 'owned_guard.py'):
    source = OWNER_RAW / name
    assert identity(source) == records[str(source.relative_to(OWNER))]
    shutil.copy2(source, WORK / name)
wrapper_path = OWNER_RAW / 'run_radiation.py'
assert identity(wrapper_path) == records[str(wrapper_path.relative_to(OWNER))]
wrapper = wrapper_path.read_text()
node = next(n for n in ast.parse(wrapper).body if isinstance(n, ast.Assign) and
            any(isinstance(t, ast.Name) and t.id == 'NAMES' for t in n.targets))
lines = wrapper.splitlines(keepends=True)
lines[node.lineno - 1:node.end_lineno] = ['NAMES = ' + repr(set(NAMES)) + '\n']
wrapper = ''.join(lines).replace('One frozen nine-case original CPU radiation action',
                               'One frozen four-case original P4 Kerr authority action')
wrapper = wrapper.replace('original-nine-cpu-radiation-controls', CASE)
assert wrapper.count('len(NAMES) == 9') == wrapper.count('len(cases) == 9') == 1
wrapper = wrapper.replace('len(NAMES) == 9', 'len(NAMES) == 4').replace('len(cases) == 9', 'len(cases) == 4')
(WORK / 'run_kerr_authority.py').write_text(wrapper)
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
p4 = next(row for row in model['acceptance_criteria'] if row['id'] == 'P4')
assert set(NAMES) <= {row['name'] for row in p4['evidence'] if row['kind'] == 'gtest'}
producers = {}
extra = []
review_path = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14/review.json'
review = json.loads(review_path.read_text())
for role, size, sha in (
    ('core', 5199120, 'c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58'),
    ('oracle', 1678696, 'be9a6f6136937567e9537e4ff25ee71589c1f98c46aa04b25790a922832abba4'),
):
    path = ROOT / ('bin/linux-gcc/tests/' + role + '/sirius_' + role + '_tests')
    recovery = ROOT / ('attestations/alignment/4bcd987/gcc-revision-synchronisation/prior-versions/' + sha)
    pin = {'bytes': size, 'sha256': sha}
    assert identity(path) == identity(recovery) == pin
    record = review['current_live_byte_bindings']['tested_artifacts']['records']['sirius_' + role + '_tests']
    assert {k: record[k] for k in ('bytes', 'sha256')} == pin
    producers[role] = {'path': key(path), 'recovery_archive': key(recovery), **pin}
    extra += [path, recovery]
snapshots = {}
for name in [*SOURCES, *HEADERS, 'tests/core/CMakeLists.txt', 'tests/oracle/CMakeLists.txt',
             'tests/operating_model.json', 'tests/labels/CTestLabels.cmake', 'tests/support/fail_on_skip.h']:
    path = ROOT / name
    assert git('show', PRODUCT + ':' + name) == path.read_bytes(), name
    snapshots[name] = {**identity(path), 'git_blob': git('rev-parse', PRODUCT + ':' + name).decode().strip()}
    target = WORK / 'reference' / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)
    pins[name], pins[key(target)] = identity(path), identity(target)
commands = json.loads((ROOT / 'bin/linux-gcc/compile_commands.json').read_text())
compiles, registrations, source_for_name = [], [], {}
for source, role in SOURCES.items():
    selected = [row for row in commands if Path(row['file']).resolve() == ROOT / source]
    assert len(selected) == 1
    row = selected[0]
    assert row['directory'] == str(ROOT / 'bin/linux-gcc')
    assert 'sirius_' + role + '_tests.dir/' in row['output'] and row['output'].endswith(Path(source).name + '.o')
    for flag in ('/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-O3', '-DNDEBUG', '-std=c++2c', '-Werror', 'tests/support/fail_on_skip.h'):
        assert flag in row['command'], (source, flag)
    compiles.append(row)
    source_names = {'{}.{}'.format(*pair) for pair in re.findall(r'TEST\(\s*(\w+)\s*,\s*(\w+)\s*\)', (ROOT / source).read_text())}
    names = sorted(set(NAMES) & source_names)
    assert len(names) == 2
    registration = ROOT / ('bin/linux-gcc/tests/' + role + '/sirius_' + role + '_tests[1]_tests.cmake')
    lines = registration.read_text().splitlines()
    for name in names:
        adds = [line for line in lines if line.startswith('add_test([=[' + name + ']=]')]
        props = [line for line in lines if line.startswith('set_tests_properties([=[' + name + ']=]')]
        assert len(adds) == len(props) == 1 and str(ROOT / producers[role]['path']) + ' ' in adds[0]
        assert '--gtest_filter=' + name + ']==]' in adds[0] and 'LABELS Correctness' in props[0] and 'TIMEOUT' not in props[0]
        assert name not in source_for_name
        source_for_name[name] = source
        registrations.append({'name': name, 'source': source, 'registration_file': key(registration),
                              'add_test': adds[0], 'properties': props[0], 'producer': producers[role]['path']})
    extra.append(registration)
assert set(source_for_name) == set(NAMES)
ninja = (ROOT / 'bin/linux-gcc/build.ninja').read_text().splitlines()
links = [line for line in ninja if line.startswith(('build tests/core/sirius_core_tests ', 'build tests/oracle/sirius_oracle_tests '))
         and ': CXX_EXECUTABLE_LINKER__' in line]
assert len(links) == 2
assert all(sum(row['output'] in line for line in links) == 1 for row in compiles)
cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
for name, value in base['source_and_registration_join']['cache_required'].items():
    assert re.search(r'^' + re.escape(name) + r':[^=]+=' + re.escape(value) + '$', cache, re.M), name
ctest = re.search(r'^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$', cache, re.M)
assert ctest
scopes = {
    NAMES[0]: {'mass_dimensionless_spin_radius_over_mass': [[0.5, -0.998, 12], [2, -0.5, 12], [3, -0.0009, 12], [0.5, 0, 12], [2, 0.0009, 12], [3, 0.5, 12], [2, 0.998, 12]],
        'assertions': 'Original optionals present and EXPECT_DOUBLE_EQ between shared authority, KerrSchildFamily ISCO, mass-scaled disk ISCO, disk angular velocity and transfer emitter angular velocity; three original extremal/zero-spin ISCO endpoint checks at mass2:12,2,18.',
        'independence_limit': 'Cross-consumer consistency against one shared CPU authority; not an independent general ISCO/orbit oracle or a check of every CPU consumer at every input.'},
    NAMES[1]: {'assertions': 'Five original invalid ISCO calls and five invalid angular-velocity calls decline; original two AccretionDiskD precondition death assertions. Inputs include nonpositive mass/radius, NaN, infinity, spin>mass and undefined retrograde denominator; no added refusal/domain boundary.'},
    NAMES[2]: {'fixture': {'mass': 1, 'spin': 0.7, 'radius': 6, 'cos_theta': 0.35, 'phi': 0.4, 'energy': 1, 'angular_momentum': 1.4, 'observer_frequency': 0.93},
        'assertions': 'Independent Boyer-Lindquist metric gives omega=-g_tphi/g_phiphi, timelike norm, u^t and u^t(E-omega*L); original four EXPECT_NEAR comparisons at2e-15. Equatorial special-case consistency against disk ZAMO uses original EXPECT_DOUBLE_EQ.',
        'independence_limit': 'One original off-equatorial physical fixture; equatorial consistency is shared-production comparison, not another independent physical oracle.'},
    NAMES[3]: {'fixture': {'mass': 1, 'spin': 0.7, 'radius': 6, 'cos_theta': 0.35, 'phi': 0.4},
        'assertions': 'Original live inverse-metric Kerr-Schild slicing normal transformed to BL has abs(radial)>1e-3 and angular rate unequal to ZAMO; transfer exists. Finite witness distinguishing radial slicing motion from a stationary azimuthal ZAMO.'},
}
assert set(scopes) == set(NAMES)
historical = ROOT / 'attestations/scientific-core/5e70607/final-build-identity-5e70607.json'
old_xml = historical.parent / 'final-integrated-cpu.xml'
history = json.loads(historical.read_text())
assert all(history['artifacts'][record['path']] != record['sha256'] for record in producers.values())
old_passes = [row.attrib for row in ET.parse(old_xml).getroot().iter('testcase') if row.get('name') in NAMES and
              not any(row.findall(tag) for tag in ('failure', 'error', 'skipped'))]
assert {row['name'] for row in old_passes} == {name for name in NAMES if name.startswith('KerrZamoTransfer.')}
extra += [Path(ctest.group(1)), Path('/usr/bin/python3.12'), review_path, OWNER / 'manifest.json',
          OWNER / 'disposition.json', OWNER_RAW / 'preparation.json', OWNER_RAW / 'run_owned.py',
          OWNER_RAW / 'owned_guard.py', wrapper_path, historical, old_xml,
          *[WORK / name for name in ('prepare.py', 'run_kerr_authority.py', 'run_owned.py', 'owned_guard.py')]]
for path in extra:
    pins[key(path)] = identity(path)
action = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'scope': 'Four original finite P4 CPU Kerr orbit/ZAMO controls only; no device/frame, all-spin/latitude/domain proof, full P4/science/Mandatory or release qualification.',
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
    'ledger_scope': {'criterion': 'P4', 'selected_original_names': NAMES, 'selected_count': 4,
        'other_completed_controls_repeated': False},
    'reference': {'original_source_snapshots': snapshots,
        'independent_reference': 'Original KerrMetricD Boyer-Lindquist metric, u_phi=0 and u.u=-1, then invariant photon contraction. No new reference implementation or threshold.',
        'scope_boundary': 'Shared CPU orbit authority consistency and explicit refusal samples plus one original arbitrary-latitude ZAMO and slicing-normal fixture; no continuous-domain mathematical proof or final integrated qualification.',
        'historical_execution': {'xml': key(old_xml), 'identity': key(historical), 'revision': history['revision'],
            'producer_sha256': {role: history['artifacts'][record['path']] for role, record in producers.items()},
            'selected_old_passes': old_passes,
            'reuse_boundary': 'Two ZAMO names pass on a different historical oracle ELF; the two orbit names are absent. Bounded review of1460 preserved XML files found no selected current-producer execution. Historical results retained as historical, not relabelled; completed radiation/Page-Thorne controls excluded.'}},
    'source_and_registration_join': {'test_sources': list(SOURCES), 'production_headers': HEADERS,
        'producers': producers,
        'accepted_source_build_evidence': base['source_and_registration_join']['accepted_source_build_evidence'],
        'cache_required': base['source_and_registration_join']['cache_required'],
        'original_compile_commands': compiles, 'original_link_ninja_lines': links,
        'registrations': registrations, 'selected_core_cases': 2, 'selected_oracle_cases': 2,
        'original_test_TIMEOUT_added': False},
    'owner_reuse': {'archive': key(OWNER), 'manifest_sha256': identity(OWNER / 'manifest.json')['sha256'],
        **{name: {'source': key(OWNER_RAW / name), **identity(WORK / name)} for name in ('run_owned.py', 'owned_guard.py')}},
    'private_tmpdir': {'path': str(WORK / 'temporary'), 'created_only_at_execution': True, 'mode': '0700',
        'owner_uid': os.getuid(), 'preserve_natural_residuals': True}, 'pins': pins})
for path in WORK.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
print(json.dumps({'prepared_without_execution': True, 'names': NAMES, 'pins': len(pins),
    'producers': producers, 'preparation': identity(WORK / 'preparation.json'),
    'wrapper': identity(WORK / 'run_kerr_authority.py'),
    'command': ['python3.12', '-B', key(WORK / 'run_kerr_authority.py')]}, indent=2))
