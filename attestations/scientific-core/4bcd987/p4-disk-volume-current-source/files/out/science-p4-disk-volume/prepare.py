"""Seal four original P4 CPU disk/volume controls without executing them."""
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
CASE = 'original-four-disk-volume-controls'
SOURCES = {'tests/core/novikov_thorne_disk_test.cpp': 'core', 'tests/render/geodesic_tracer_test.cpp': 'render'}
NAMES = ['AccretionDiskTest.ShakuraSunyaevProfileHasZeroTorqueEdgeAndDeclaredScale', 'GeodesicTracerVolumetric.OpticallyThinTransferIsNotDiscardedAtCompositionBoundary', 'GeodesicTracerVolumetric.RedshiftAndDopplerReachTheLiveVolumeSource', 'GeodesicTracerVolumetric.TransferAccumulatesAcrossEveryTraversedSegment']
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
                               'One frozen four-case original P4 disk/volume action')
wrapper = wrapper.replace('original-nine-cpu-radiation-controls', CASE)
assert wrapper.count('len(NAMES) == 9') == wrapper.count('len(cases) == 9') == 1
wrapper = wrapper.replace('len(NAMES) == 9', 'len(NAMES) == 4').replace('len(cases) == 9', 'len(cases) == 4')
(WORK / 'run_disk_volume.py').write_text(wrapper)
model = json.loads((ROOT / 'tests/operating_model.json').read_text())
p4 = next(row for row in model['acceptance_criteria'] if row['id'] == 'P4')
assert set(NAMES) <= {row['name'] for row in p4['evidence'] if row['kind'] == 'gtest'}
producers = {}
extra = []
review_path = ROOT / 'attestations/fresh-session-review/4bcd987/2026-10-14/review.json'
review = json.loads(review_path.read_text())
for role, size, sha in (
    ('core', 5199120, 'c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58'),
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
for name in [*SOURCES, *HEADERS, 'tests/core/CMakeLists.txt', 'tests/render/CMakeLists.txt',
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
    assert len(names) == {'core': 1, 'render': 3}[role]
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
links = [line for line in ninja if line.startswith(('build tests/core/sirius_core_tests ', 'build src/sirius/app/sirius_render_tests '))
         and ': CXX_EXECUTABLE_LINKER__' in line]
assert len(links) == 2
assert all(sum(row['output'] in line for line in links) == 1 for row in compiles)
cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
for name, value in base['source_and_registration_join']['cache_required'].items():
    assert re.search(r'^' + re.escape(name) + r':[^=]+=' + re.escape(value) + '$', cache, re.M), name
ctest = re.search(r'^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$', cache, re.M)
assert ctest
scopes = {'AccretionDiskTest.ShakuraSunyaevProfileHasZeroTorqueEdgeAndDeclaredScale': {'fixture': {'inner_radius': 6, 'temperature_scale': 10000, 'reference_radius': 9, 'sample_radius': 36}, 'assertions': 'Original EXPECT_DOUBLE_EQ: edge flux/temperature zero, below-edge temperature zero, scale at1.5*r_in, NaN scale returns zero. Original independently expanded Newtonian zero-torque normalized temperature expression at r36 within1e-12 and below scale.', 'independence_limit': 'One declared radial profile fixture; no relativistic or alpha-disk vertical-structure claim.'}, 'GeodesicTracerVolumetric.OpticallyThinTransferIsNotDiscardedAtCompositionBoundary': {'fixture': {'escape_radius': 80, 'disk_temperature_inner': 1, 'disk_temperature_model': 'ShakuraSunyaev', 'disk_outer': 20, 'camera_size': [3, 3], 'pixel': [1, 1], 'sample': [0.5, 0.5], 'camera_theta': 'pi/2', 'max_step': 1, 'mass': 1, 'spin': 0, 'max_steps': 2000, 'disk_inner': 6, 'height_ratio': 0.1, 'tau_midplane': 1e-05, 'volume_samples': 4, 'initial_step': 0.5, 'camera_radius': 50}, 'rays': 1, 'assertions': 'Original positive optical depth strictly below0.01, volumetric_hit true and first emission channel positive. No additional outcome or numerical_failure assertion.'}, 'GeodesicTracerVolumetric.RedshiftAndDopplerReachTheLiveVolumeSource': {'fixture': {'escape_radius': 80, 'disk_temperature_inner': 1, 'disk_temperature_model': 'ShakuraSunyaev', 'disk_outer': 20, 'camera_size': [3, 3], 'pixel': [1, 1], 'sample': [0.5, 0.5], 'camera_theta': 'pi/2', 'max_step': 1, 'mass': 1, 'spin': 0.6, 'max_steps': 4000, 'disk_inner': 4, 'height_ratio': 0.15, 'tau_midplane': 1, 'volume_samples': 2, 'initial_step': 0.25, 'camera_radius': 40}, 'rays': 2, 'assertions': 'Same original ray with doppler_beaming true/false: both volume hits and positive first emission channels; beamed first channel strictly below diagnostic gravity-only/ZAMO branch. Finite live-consumer ordering, not an independent absolute radiation oracle.'}, 'GeodesicTracerVolumetric.TransferAccumulatesAcrossEveryTraversedSegment': {'fixture': {'escape_radius': 80, 'disk_temperature_inner': 1, 'disk_temperature_model': 'ShakuraSunyaev', 'disk_outer': 20, 'camera_size': [3, 3], 'pixel': [1, 1], 'sample': [0.5, 0.5], 'camera_theta': 'pi/2', 'max_step': 1, 'mass': 1, 'spin': 0, 'max_steps': 20000, 'horizon_factor': 1.05, 'disk_inner': 6, 'height_ratio': 0.1, 'tau_midplane': 2, 'volume_samples': 4, 'initial_step': 0.5, 'camera_radius': 50}, 'rays': 2, 'assertions': 'Original uncapped outcome Horizon, volume hit, active affine length>10, optical depth>0.1, finite positive first emission channel. Original second ray tau ceiling0.25: volume hit, optical depth within1e-6 and <=0.25+1e-6. No new pointwise/segmentwise oracle.'}}
assert set(scopes) == set(NAMES)
reuse_root = ROOT / 'attestations/software-vulkan/ef88b99/retained-dopri-integration'
reuse = reuse_root / 'diagnostic'
reuse_manifest = json.loads((reuse_root / 'manifest.json').read_text())
reuse_records = {row['path']: {k: row[k] for k in ('bytes', 'sha256')} for row in reuse_manifest['files']}
reuse_paths = [reuse / name for name in ('host-controls/report.json', 'host-controls-inputs-before.json',
    'host-controls-inputs-after.json', 'host-controls-junit.xml', 'host-controls-LastTest.log')]
for path in reuse_paths:
    assert identity(path) == reuse_records[str(path.relative_to(reuse_root))]
reuse_before = json.loads((reuse / 'host-controls-inputs-before.json').read_text())
assert reuse_before == json.loads((reuse / 'host-controls-inputs-after.json').read_text())
reuse_core = next(row for row in reuse_before if row['resolved_path'] == str(ROOT / producers['core']['path']))
assert {k: reuse_core[k] for k in ('bytes', 'sha256')} == {k: producers['core'][k] for k in ('bytes', 'sha256')}
reused_name = 'VolumetricDiskClosure.TruncatedGaussianColumnEqualsDeclaredOpticalDepth'
reuse_case = next(row for row in ET.parse(reuse / 'host-controls-junit.xml').getroot().iter('testcase') if row.get('name') == reused_name)
assert reuse_case.get('status') == 'run' and not any(reuse_case.findall(tag) for tag in ('failure', 'error', 'skipped'))
reuse_owner = json.loads((reuse / 'host-controls/report.json').read_text())
assert reuse_owner['status'] == 'terminal' and reuse_owner['passed'] and reuse_owner['child_exit'] == 0
assert reuse_owner['source_unchanged'] and reuse_owner['stop_reason'] is None and not reuse_owner['remaining_owned_processes']
for name in ('tests/core/volumetric_disk_test.cpp', 'src/sirius/core/disk/volumetric_disk.h'):
    assert git('show', 'ef88b99b74b6741ecfb6e009403d9ea10cef26e4:' + name) == (ROOT / name).read_bytes()
production_compile = [row for row in commands if Path(row['file']).resolve() == ROOT / 'src/sirius/backend/cpu/geodesic_tracer.cpp']
assert len(production_compile) == 1
for flag in ('/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-O3', '-DNDEBUG', '-std=c++2c', '-Werror'):
    assert flag in production_compile[0]['command']
production_archive = [line for line in ninja if line.startswith('build src/sirius/backend/libsirius_backend_cpu.a:')]
assert len(production_archive) == 1 and production_compile[0]['output'] in production_archive[0]
assert any(line.startswith('build src/sirius/app/sirius_render_tests ') and 'src/sirius/backend/libsirius_backend_cpu.a' in line for line in links)
extra += [Path(ctest.group(1)), Path('/usr/bin/python3.12'), review_path, OWNER / 'manifest.json',
          OWNER / 'disposition.json', OWNER_RAW / 'preparation.json', OWNER_RAW / 'run_owned.py',
          OWNER_RAW / 'owned_guard.py', wrapper_path, reuse_root / 'manifest.json', *reuse_paths,
          *[WORK / name for name in ('prepare.py', 'run_disk_volume.py', 'run_owned.py', 'owned_guard.py')]]
for path in extra:
    pins[key(path)] = identity(path)
action = {'revision': RUNNER, 'kind': 'check', 'expected_exit': 0,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'scope': 'Four original finite P4 CPU Shakura-Sunyaev/volume-transfer controls only; five explicit CPU rays, no device/frame, exhaustive transfer proof, full P4/science/Mandatory or release qualification.',
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
        'independent_reference': 'Original expanded Newtonian zero-torque flux/temperature expression. Three original live CPU volume tests provide finite accumulation/cap, thin-layer retention and Doppler-branch ordering assertions; no new implementation or threshold.',
        'cost_scope': {'full_frames': 0, 'explicit_cpu_rays': 5, 'sum_original_trace_iteration_caps': 50000,
            'bounded_sampling': 'Original volume samples4/4/2 per accepted segment; rejected trace attempts also count against max_steps. No wall-time forecast or completion guarantee; outer diagnostic cap300s unchanged.'},
        'scope_boundary': 'Finite phenomenological radial/volume closure and live-consumer witnesses; not an independent absolute volume-radiation or pointwise proof.',
        'matching_historical_reuse': {'name': reused_name, 'xml': key(reuse / 'host-controls-junit.xml'),
            'owner': key(reuse / 'host-controls/report.json'), 'input_before': key(reuse / 'host-controls-inputs-before.json'),
            'input_after': key(reuse / 'host-controls-inputs-after.json'), 'manifest': key(reuse_root / 'manifest.json'),
            'case': reuse_case.attrib, 'same_current_core_producer': reuse_core,
            'source_join': 'Original test and closure header byte-equal between ef88b99 and current product4bcd987; producer matches accepted current core exactly.',
            'reuse_scope': 'Reuse the completed finite Gaussian-column assertion; historical owner used1200s/12GiB/2GiB and is not relabelled as the new300s/4GiB action.'},
        'gap_boundary': 'Bounded1461-preserved-XML review found no selected four-case execution; excludes passed radiation/Page-Thorne/Kerr controls and same-producer Gaussian-column reuse.'},
    'source_and_registration_join': {'test_sources': list(SOURCES), 'production_headers': HEADERS,
        'producers': producers,
        'accepted_source_build_evidence': base['source_and_registration_join']['accepted_source_build_evidence'],
        'cache_required': base['source_and_registration_join']['cache_required'],
        'original_compile_commands': compiles, 'original_link_ninja_lines': links,
        'registrations': registrations, 'selected_core_cases': 1, 'selected_render_cases': 3,
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
    'wrapper': identity(WORK / 'run_disk_volume.py'),
    'command': ['python3.12', '-B', key(WORK / 'run_disk_volume.py')]}, indent=2))
