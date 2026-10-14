"""Seal four original remaining P4 CPU non-frame controls without executing them."""
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
CASE = 'original-four-non-frame-cpu-controls'
SOURCES = {'tests/backend/cpu_trace_boundary_test.cpp': 'backend', 'tests/render/geodesic_tracer_test.cpp': 'render'}
NAMES = ['CpuTraceBoundary.TruncatedPageThorneLiveProfileUsesDeclaredZeroTorqueEdge', 'GeodesicTracerTest.LiveDiskTemperatureUsesFullPageThorneProfile', 'GeodesicTracerTest.LiveDiskTemperatureUsesZeroTorqueShakuraSunyaevProfile', 'GeodesicTracerRedshift.NearExtremalInnerDiskEmissionRemainsFinite']
HEADERS = ['src/sirius/core/disk/novikov_thorne_disk.h', 'src/sirius/core/disk/disk_model.h', 'src/sirius/core/disk/volumetric_disk.h', 'src/sirius/core/constants.h', 'src/sirius/core/relativistic_transfer.h', 'src/sirius/core/kerr_orbits.h', 'src/sirius/core/coordinates.h', 'src/sirius/core/metrics/kerr_schild_family.h', 'src/sirius/core/metrics/outgoing_kerr_schild.h', 'src/sirius/core/camera.h', 'src/sirius/core/camera_launch.h', 'src/sirius/core/observer_frame.h', 'src/sirius/core/geodesic_integrator.h', 'src/sirius/core/spectral/colour_modes.h', 'src/sirius/core/spectral/spectral_radiance.h', 'src/sirius/core/disk/turbulence.h', 'src/sirius/backend/cpu/geodesic_tracer.h', 'src/sirius/backend/cpu/geodesic_tracer.cpp', 'src/sirius/backend/trace_step_executor.h', 'src/sirius/backend/CMakeLists.txt', 'src/sirius/render/trace_domain.h', 'tests/render/doppler_toggle_test.cpp', 'src/sirius/base/contracts.h']
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
                               'One frozen four-case original P4 non-frame remainder action')
wrapper = wrapper.replace('original-nine-cpu-radiation-controls', CASE)
assert wrapper.count('len(NAMES) == 9') == wrapper.count('len(cases) == 9') == 1
wrapper = wrapper.replace('len(NAMES) == 9', 'len(NAMES) == 4').replace('len(cases) == 9', 'len(cases) == 4')
(WORK / 'run_remainder.py').write_text(wrapper)
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
p4 = next(row for row in model['acceptance_criteria'] if row['id'] == 'P4')
assert set(NAMES) <= {row['name'] for row in p4['evidence'] if row['kind'] == 'gtest'}
producers = {}
extra = []
review_path = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14/review.json'
review = json.loads(review_path.read_text())
for role, size, sha in (
    ('backend', 13440952, 'e5f3b0f120ec6a12490e3e9b6e491d19eb12bfad62c1cdce11722c091a152149'),
    ('render', 12208240, '2489bc56bd0c1c53ab01da021119d6292d8fdb6f51037b82819d3dcdf6e87ab4'),
):
    path = ROOT / ('bin/linux-gcc/src/sirius/app/sirius_render_tests' if role == 'render' else 'bin/linux-gcc/tests/backend/sirius_backend_tests')
    recovery = ROOT / ('attestations/alignment/4bcd987/gcc-revision-synchronisation/' + ('current-versions/') + sha)
    pin = {'bytes': size, 'sha256': sha}
    assert identity(path) == identity(recovery) == pin
    record = review['current_live_byte_bindings']['tested_artifacts']['records']['sirius_' + role + '_tests']
    assert {k: record[k] for k in ('bytes', 'sha256')} == pin
    producers[role] = {'path': key(path), 'recovery_archive': key(recovery), **pin}
    extra += [path, recovery]
snapshots = {}
for name in [*SOURCES, *HEADERS, 'tests/backend/CMakeLists.txt', 'tests/render/CMakeLists.txt',
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
    source_names = {'{}.{}'.format(*pair) for pair in re.findall(r'TEST(?:_F)?\(\s*(\w+)\s*,\s*(\w+)\s*\)', (ROOT / source).read_text())}
    names = sorted(set(NAMES) & source_names)
    assert len(names) == {'backend': 1, 'render': 3}[role]
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
links = [line for line in ninja if line.startswith(('build tests/backend/sirius_backend_tests ', 'build src/sirius/app/sirius_render_tests '))
         and ': CXX_EXECUTABLE_LINKER__' in line]
assert len(links) == 2
assert all(sum(row['output'] in line for line in links) == 1 for row in compiles)
cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
for name, value in base['source_and_registration_join']['cache_required'].items():
    assert re.search(r'^' + re.escape(name) + r':[^=]+=' + re.escape(value) + '$', cache, re.M), name
ctest = re.search(r'^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$', cache, re.M)
assert ctest
scopes = {'CpuTraceBoundary.TruncatedPageThorneLiveProfileUsesDeclaredZeroTorqueEdge': {'original_source': 'tests/backend/cpu_trace_boundary_test.cpp:601', 'fixture': {'mass': 1, 'spin': 0, 'camera_size': [64, 64], 'camera_radius': 50, 'camera_theta': 'pi/2', 'pixel_range': 'x,y=8,10,...,54', 'sample': [0.5, 0.5], 'disk_outer': 20, 'old_inner_edge': 6, 'new_inner_edge': 10, 'temperature_scale': 10000, 'max_steps_per_trace': 5000}, 'search': 'One original cache-populating ray(32,32), SetConfig changes the declared inner edge to10, then original raster search until compared>=12 or exhausted; at most577 Trace calls.', 'assertions': 'Original reference temperature positive; each valid crossing temperature agrees within2e-3 with scale*AccretionDiskD.Temperature(r)/Temperature(1.5*InnerRadius), and compared>=6.', 'reference_limit': 'Live consumer/cache invalidation check against the production AccretionDiskD authority at a declared truncated edge, not an independent Page-Thorne quadrature.'}, 'GeodesicTracerTest.LiveDiskTemperatureUsesFullPageThorneProfile': {'original_source': 'tests/render/geodesic_tracer_test.cpp:364', 'fixture': {'mass': 1, 'spin': 0, 'camera_size': [64, 64], 'camera_radius': 50, 'camera_theta': 'pi/2', 'pixel_range': 'x,y=8,10,...,54', 'sample': [0.5, 0.5], 'disk_inner': 6, 'disk_outer': 20, 'max_steps_per_trace': 5000}, 'search': 'Unchanged fixture raster search until compared>=12 or exhausted; at most576 Trace calls.', 'assertions': 'Original positive reference temperature; each valid crossing temperature agrees within2e-5 with AccretionDiskD.Temperature(r)/Temperature(1.5*IscoRadius), and compared>=6.', 'reference_limit': 'Live default-profile consumer against shared production AccretionDiskD authority; independent flux quadrature is covered separately by the completed Page-Thorne controls.'}, 'GeodesicTracerTest.LiveDiskTemperatureUsesZeroTorqueShakuraSunyaevProfile': {'original_source': 'tests/render/geodesic_tracer_test.cpp:338', 'fixture': {'mass': 1, 'spin': 0, 'camera_size': [64, 64], 'camera_radius': 50, 'camera_theta': 'pi/2', 'pixel_range': 'x,y=8,10,...,54', 'sample': [0.5, 0.5], 'disk_inner': 6, 'disk_outer': 20, 'temperature_scale': 10000, 'max_steps_per_trace': 5000}, 'search': 'Unchanged fixture raster search until compared>=12 or exhausted; at most576 Trace calls.', 'assertions': 'Original each valid crossing temperature agrees within2e-3 with shared ShakuraSunyaevTemperature(scale,r,inner), and compared>=6.', 'reference_limit': 'Live selected-profile consumer consistency against the production radial helper; no independent radiation or vertical-structure oracle.'}, 'GeodesicTracerRedshift.NearExtremalInnerDiskEmissionRemainsFinite': {'original_source': 'tests/render/geodesic_tracer_test.cpp:727', 'fixture': {'mass': 1, 'spin': 0.998, 'camera_size': [96, 54], 'camera_radius': 30, 'camera_theta': '80*pi/180', 'fov_degrees': 55, 'pixel_range': 'x=4,6,...,90;y=4,6,...,48', 'sample': [0.5, 0.5], 'disk_inner': 'AccretionDiskD.ComputeIsco(0.998)', 'disk_outer': 20, 'initial_step': 0.1, 'max_step': 1, 'max_steps_per_trace': 20000}, 'search': 'Original raster search until inner_crossings>=3 or exhausted; at most1012 Trace calls.', 'assertions': 'For each valid crossing with r<2M, original redshift finite and positive; inner_crossings>=1. No extra termination or numerical_failure assertion added.', 'reference_limit': 'Finite inner-disk emission witness; no independent absolute transfer value or all-ray proof.'}}
assert set(scopes) == set(NAMES)
previous = ROOT / 'attestations/scientific-core/4bcd987/p4-disk-volume-current-source'
previous_raw = previous / 'files/out/science-p4-disk-volume'
previous_records = json.loads((previous / 'manifest.json').read_text())['files']
for name in ('prepare.py', 'preparation.json'):
    relative = str((previous_raw / name).relative_to(previous))
    assert identity(previous_raw / name) == previous_records[relative]
previous_preparation = json.loads((previous_raw / 'preparation.json').read_text())
reused_closure = previous_preparation['reference']['matching_historical_reuse']
production_compile = [row for row in commands if Path(row['file']).resolve() == ROOT / 'src/sirius/backend/cpu/geodesic_tracer.cpp']
assert len(production_compile) == 1
for flag in ('/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-O3', '-DNDEBUG', '-std=c++2c', '-Werror'):
    assert flag in production_compile[0]['command']
production_archive = [line for line in ninja if line.startswith('build src/sirius/backend/libsirius_backend_cpu.a:')]
assert len(production_archive) == 1 and production_compile[0]['output'] in production_archive[0]
assert any(line.startswith('build src/sirius/app/sirius_render_tests ') and 'src/sirius/backend/libsirius_backend_cpu.a' in line for line in links)
extra += [Path(ctest.group(1)), Path('/usr/bin/python3.12'), review_path, OWNER / 'manifest.json',
          OWNER / 'disposition.json', OWNER_RAW / 'preparation.json', OWNER_RAW / 'run_owned.py',
          OWNER_RAW / 'owned_guard.py', wrapper_path, previous / 'manifest.json', previous_raw / 'prepare.py', previous_raw / 'preparation.json',
          *[WORK / name for name in ('prepare.py', 'run_remainder.py', 'run_owned.py', 'owned_guard.py')]]
for path in extra:
    pins[key(path)] = identity(path)
action = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'scope': 'Four unchanged original finite P4 CPU live-temperature/inner-disk controls only; original central-ray search loops and thresholds, no frame publication/device execution, full P4/science/Mandatory or release qualification.',
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
        'independent_reference': 'Original live-temperature profile/cache checks consume shared production AccretionDiskD/ShakuraSunyaev authorities; the near-extremal scan asserts finite positive inner-disk redshift. Independent flux/radiation controls remain separate completed evidence; no new implementation or threshold.',
        'cost_scope': {'full_frame_publications': 0, 'maximum_original_Trace_calls': 2741, 'maximum_sum_trace_attempt_caps': 28885000, 'scope': 'These are source loop/attempt upper bounds, not timing predictions. Original searches stop by their existing crossing counters. Outer300s/4GiB/2GiB limits unchanged; timeout means incomplete.'},
        'scope_boundary': 'Finite live profile-consumer/cache and inner-disk redshift assertions only; no independent general radiation proof.',
        'matching_historical_reuse': reused_closure,
        'gap_boundary': 'Bounded preserved-XML name review found no selected execution on matching current producers. Only CpuTraceBoundary has three historical XML matches on older backend producers; none is relabelled. The original two-pass120x68 DopplerToggleTest.SuppressionCollapsesDiskAsymmetry remains pending as a separate later owned action; no one-ray substitution.'},
    'source_and_registration_join': {'test_sources': list(SOURCES), 'production_headers': HEADERS,
        'producers': producers,
        'accepted_source_build_evidence': base['source_and_registration_join']['accepted_source_build_evidence'],
        'cache_required': base['source_and_registration_join']['cache_required'],
        'original_compile_commands': compiles, 'original_link_ninja_lines': links,
        'registrations': registrations, 'selected_backend_cases': 1, 'selected_render_cases': 3,
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
    'wrapper': identity(WORK / 'run_remainder.py'),
    'command': ['python3.12', '-B', key(WORK / 'run_remainder.py')]}, indent=2))
