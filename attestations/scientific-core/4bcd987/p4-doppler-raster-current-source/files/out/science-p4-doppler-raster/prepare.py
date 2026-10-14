"""Seal the original two-pass P4 CPU Doppler raster without executing it."""
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
CASE = 'original-doppler-raster-control'
SOURCES = {'tests/render/doppler_toggle_test.cpp': 'render'}
NAMES = ['DopplerToggleTest.SuppressionCollapsesDiskAsymmetry']
HEADERS = ['src/sirius/core/disk/novikov_thorne_disk.h', 'src/sirius/core/disk/disk_model.h', 'src/sirius/core/disk/volumetric_disk.h', 'src/sirius/core/constants.h', 'src/sirius/core/relativistic_transfer.h', 'src/sirius/core/kerr_orbits.h', 'src/sirius/core/coordinates.h', 'src/sirius/core/metrics/kerr_schild_family.h', 'src/sirius/core/metrics/outgoing_kerr_schild.h', 'src/sirius/core/camera.h', 'src/sirius/core/camera_launch.h', 'src/sirius/core/observer_frame.h', 'src/sirius/core/geodesic_integrator.h', 'src/sirius/core/spectral/colour_modes.h', 'src/sirius/core/spectral/spectral_radiance.h', 'src/sirius/core/disk/turbulence.h', 'src/sirius/backend/cpu/geodesic_tracer.h', 'src/sirius/backend/cpu/geodesic_tracer.cpp', 'src/sirius/backend/trace_step_executor.h', 'src/sirius/backend/CMakeLists.txt', 'src/sirius/base/contracts.h']
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
adapted_preparer = ROOT / 'attestations/scientific-core/4bcd987/p4-disk-volume-current-source/files/out/science-p4-disk-volume/prepare.py'
assert identity(adapted_preparer)['sha256'] == 'e2e36d2c85468956b4eb0936920d75dd54995185086fde5d74c6f7f1b6562b49'
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
assert wrapper.count('LIMITS = {\"seconds\": 300,') == 1
wrapper = wrapper.replace('LIMITS = {\"seconds\": 300,', 'LIMITS = {\"seconds\": 900,', 1)
node = next(n for n in ast.parse(wrapper).body if isinstance(n, ast.Assign) and
            any(isinstance(t, ast.Name) and t.id == 'NAMES' for t in n.targets))
lines = wrapper.splitlines(keepends=True)
lines[node.lineno - 1:node.end_lineno] = ['NAMES = ' + repr(set(NAMES)) + '\n']
wrapper = ''.join(lines).replace('One frozen nine-case original CPU radiation action',
                               'One frozen one-case original P4 Doppler raster action')
wrapper = wrapper.replace('original-nine-cpu-radiation-controls', CASE)
assert wrapper.count('len(NAMES) == 9') == wrapper.count('len(cases) == 9') == 1
wrapper = wrapper.replace('len(NAMES) == 9', 'len(NAMES) == 1').replace('len(cases) == 9', 'len(cases) == 1')
(WORK / 'run_doppler.py').write_text(wrapper)
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
p4 = next(row for row in model['acceptance_criteria'] if row['id'] == 'P4')
assert set(NAMES) <= {row['name'] for row in p4['evidence'] if row['kind'] == 'gtest'}
producers = {}
extra = []
review_path = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14/review.json'
review = json.loads(review_path.read_text())
for role, size, sha in (
    ('render', 12208240, '2489bc56bd0c1c53ab01da021119d6292d8fdb6f51037b82819d3dcdf6e87ab4'),
):
    path = ROOT / ('bin/linux-gcc/src/sirius/app/sirius_render_tests' if role == 'render' else 'bin/linux-gcc/tests/core/sirius_core_tests')
    recovery = ROOT / ('attestations/alignment/4bcd987/gcc-revision-synchronisation/' + ('current-versions/' if role == 'render' else 'prior-versions/') + sha)
    pin = {'bytes': size, 'sha256': sha}
    assert identity(path) == identity(recovery) == pin
    record = review['current_live_byte_bindings']['tested_artifacts']['records']['sirius_' + role + '_tests']
    assert {k: record[k] for k in ('bytes', 'sha256')} == pin
    producers[role] = {'path': key(path), 'recovery_archive': key(recovery), **pin}
    extra += [path, recovery]
snapshots = {}
for name in [*SOURCES, *HEADERS, 'tests/render/CMakeLists.txt',
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
    assert len(names) == {'render': 1}[role]
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
links = [line for line in ninja if line.startswith('build src/sirius/app/sirius_render_tests ')
         and ': CXX_EXECUTABLE_LINKER__' in line]
assert len(links) == 1
assert all(sum(row['output'] in line for line in links) == 1 for row in compiles)
cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
for name, value in base['source_and_registration_join']['cache_required'].items():
    assert re.search(r'^' + re.escape(name) + r':[^=]+=' + re.escape(value) + '$', cache, re.M), name
ctest = re.search(r'^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$', cache, re.M)
assert ctest
scopes = {'DopplerToggleTest.SuppressionCollapsesDiskAsymmetry': {'original_source': 'tests/render/doppler_toggle_test.cpp:134', 'fixture': {'mass': 1, 'spin': 0.9, 'camera_width': 120, 'camera_height': 68, 'passes': ['doppler_beaming=true', 'doppler_beaming=false'], 'rays_per_pass': 8160, 'total_rays': 16320, 'pixel_sample': [0.5, 0.5], 'camera_radius': 30, 'camera_theta': '80*pi/180', 'camera_phi': 0, 'fov_degrees': 55, 'disk_inner': 'AccretionDiskD.ComputeIsco(0.9)', 'disk_outer': 20, 'disk_temperature_inner': 1, 'escape_radius': 100, 'horizon_factor': 1.05, 'max_steps_per_trace': 20000, 'initial_step': 0.1, 'max_step': 2}, 'parallelism': 'Original std::clamp(std::thread::hardware_concurrency(),1u,16u); each async worker owns metric/tracer/camera, every pixel retained by raster index, all futures joined before reduction. No worker-count override.', 'image_join': 'Two complete120x68 rasters. For DiskHit pixels, reduce g^4*T^4 luminosity and (redshift/zamo_disk_redshift)^4 per-hit factors separately by image half, in original raster order. Normalize total luminosity asymmetry and mean per-hit factor asymmetry using original guards.', 'assertions': ['ASSERT_GT(on.disk_hits,100)', 'EXPECT_EQ(on.disk_hits,off.disk_hits)', 'EXPECT_GT(on.doppler_factor,0.1)', 'EXPECT_LT(off.doppler_factor,1e-6)', 'EXPECT_LT(off.observed,0.25*on.observed)'], 'reference_limit': 'Original finite live CPU raster/shading-toggle witness; the suppressed branch is the diagnostic ZAMO substitution. Not an independent absolute image-radiance oracle, all-scene bound or standalone platform qualification.', 'console_observations': 'Original [doppler] row preserves observed_on,observed_off,factor_on,factor_off,hits in verbose CTest log and JUnit system-out. Original GTest assertions remain the numerical authority.'}}
assert set(scopes) == set(NAMES)
production_compile = [row for row in commands if Path(row['file']).resolve() == ROOT / 'src/sirius/backend/cpu/geodesic_tracer.cpp']
assert len(production_compile) == 1
for flag in ('/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-O3', '-DNDEBUG', '-std=c++2c', '-Werror'):
    assert flag in production_compile[0]['command']
production_archive = [line for line in ninja if line.startswith('build src/sirius/backend/libsirius_backend_cpu.a:')]
assert len(production_archive) == 1 and production_compile[0]['output'] in production_archive[0]
assert any(line.startswith('build src/sirius/app/sirius_render_tests ') and 'src/sirius/backend/libsirius_backend_cpu.a' in line for line in links)
extra += [adapted_preparer, Path(ctest.group(1)), Path('/usr/bin/python3.12'), review_path, OWNER / 'manifest.json',
          OWNER / 'disposition.json', OWNER_RAW / 'preparation.json', OWNER_RAW / 'run_owned.py',
          OWNER_RAW / 'owned_guard.py', wrapper_path,
          *[WORK / name for name in ('prepare.py', 'run_doppler.py', 'run_owned.py', 'owned_guard.py')]]
for path in extra:
    pins[key(path)] = identity(path)
action = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0,
    'limits': {'seconds': 900, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'scope': 'One original finite P4 CPU Doppler witness with both complete120x68 rasters, original worker cap, deterministic raster reductions and all original numerical assertions. No dimensions, rays or thresholds changed; full P4/science/Mandatory or release qualification remains separate.',
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
    'adapted_preparation_source': {'path': key(adapted_preparer), **identity(adapted_preparer)},
    'ledger_scope': {'criterion': 'P4', 'selected_original_names': NAMES, 'selected_count': 1,
        'other_completed_controls_repeated': False},
    'reference': {'original_source_snapshots': snapshots,
        'independent_reference': 'Original two full camera rasters compare beaming-on circular-emitter transfer with the Doppler-suppressed diagnostic ZAMO branch, using original g^4*T^4 image asymmetry and mean isolated-factor assertions. No new reference implementation or threshold.',
        'cost_scope': {'original_full_raster_passes': 2, 'width': 120, 'height': 68, 'Trace_calls': 16320, 'maximum_sum_trace_attempt_caps': 326400000, 'workers': 'Original hardware_concurrency clamped1..16', 'original_test_TIMEOUT': None, 'scope': 'Original registration has no TIMEOUT;900s is a finite outer diagnostic cap, not an original scientific bound or a completion forecast. Preserve full work and all original assertions; timeout means incomplete.'},
        'scope_boundary': 'Finite live CPU raster Doppler-suppression assertion only; not an independent absolute image-radiance or general scene proof.',
        'gap_boundary': 'Bounded preserved-XML review found no execution of this exact original name on a matching current producer. Previously passed radiation/Page-Thorne/Kerr/disk-volume/Gaussian-closure controls are not selected or repeated.'},
    'source_and_registration_join': {'test_sources': list(SOURCES), 'production_headers': HEADERS,
        'producers': producers,
        'accepted_source_build_evidence': base['source_and_registration_join']['accepted_source_build_evidence'],
        'cache_required': base['source_and_registration_join']['cache_required'],
        'original_compile_commands': compiles, 'original_link_ninja_lines': links,
        'registrations': registrations, 'selected_render_cases': 1,
        'live_cpu_production_compile': production_compile, 'live_cpu_archive_link': production_archive,
        'original_test_TIMEOUT_added': False},
    'owner_reuse': {'archive': key(OWNER), 'manifest_sha256': identity(OWNER / 'manifest.json')['sha256'],
        **{name: {'source': key(OWNER_RAW / name), **identity(WORK / name)} for name in ('run_owned.py', 'owned_guard.py')}},
    'private_tmpdir': {'path': str(WORK / 'temporary'), 'created_only_at_execution': True, 'mode': '0700',
        'owner_uid': os.getuid(), 'preserve_natural_residuals': True}, 'pins': pins})
for path in WORK.glob('*.py'):
    ast.parse(path.read_text(), filename=str(path))
print(json.dumps({'prepared_without_execution': True, 'names': NAMES, 'pins': len(pins),
    'producers': producers, 'preparation': identity(WORK / 'preparation.json'),
    'wrapper': identity(WORK / 'run_doppler.py'),
    'command': ['python3.12', '-B', key(WORK / 'run_doppler.py')]}, indent=2))
