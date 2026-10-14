"""One plain historical export and bounded input sealing; never builds or renders."""
from pathlib import Path, PurePosixPath
import datetime
import hashlib
import io
import json
import shutil
import subprocess
import sys
import tarfile

sys.dont_write_bytecode = True
work = Path(__file__).resolve().parent
root = work.parent.parent
assert root == Path('/home/astra/.project/Sirius') and Path.cwd() == root
assert not (work / 'source').exists(), 'plain export is one-shot; never replace it silently'
source_revision = 'b739f920c2ae67847ed7961e011d2350a839e410'
legacy_revision = 'ec2d8e506dfd0477e3c9eb4b2b4aa1bee9831fe2'
original_revision = '1cc9e2a4e44de5276ae6ed5403beca7550490e85'
build = root / 'bin/linux-gcc/forensic-reconstruction-b739f92'
source = work / 'source'
reference = work / 'reference'
reference.mkdir()

def query(argv):
    return subprocess.check_output(argv, cwd=root, timeout=60)

def git(*args):
    return query(['git', *args])

def write(path, value):
    assert path.is_relative_to(work) and not path.exists(), path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

def identity(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def seal(records):
    return hashlib.sha256(json.dumps(records, sort_keys=True, separators=(',', ':')).encode()).hexdigest()

entries = {}
for item in git('ls-tree', '-rz', '--full-tree', source_revision).split(b'\0'):
    if not item:
        continue
    metadata, name = item.split(b'\t', 1)
    mode, kind, blob = metadata.decode().split()
    name = name.decode()
    assert mode == '100644' and kind == 'blob'
    assert not PurePosixPath(name).is_absolute() and '..' not in PurePosixPath(name).parts and '.git' not in PurePosixPath(name).parts
    entries[name] = {'git_blob': blob, 'git_mode': mode}
assert len(entries) == 1162
archive = git('archive', '--format=tar', source_revision)
archive_sha = hashlib.sha256(archive).hexdigest()
source.mkdir()
with tarfile.open(fileobj=io.BytesIO(archive), mode='r:') as stream:
    actual = set()
    for member in stream.getmembers():
        name = member.name.rstrip('/')
        assert not PurePosixPath(name).is_absolute() and '..' not in PurePosixPath(name).parts and '.git' not in PurePosixPath(name).parts
        if member.isdir():
            (source / name).mkdir(parents=True, exist_ok=True)
            continue
        assert member.isfile() and name in entries and name not in actual and member.size <= 31_000_000
        data = stream.extractfile(member).read()
        assert len(data) == member.size
        assert hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest() == entries[name]['git_blob']
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(0o644)
        entries[name].update(identity(path))
        actual.add(name)
assert actual == set(entries) and sum(e['bytes'] for e in entries.values()) == 43_614_536
write(work / 'source-manifest.json', {'revision': source_revision, 'tree': git('rev-parse', source_revision + '^{tree}').decode().strip(),
      'git_archive_tar_sha256': archive_sha, 'files': entries, 'file_count': len(entries), 'bytes': 43_614_536,
      'file_records_sha256': seal(entries), 'plain_export_without_git': True})

legacy_paths = ['src/Sirius.Core', 'src/Sirius.Render', 'src/Sirius.Infrastructure', 'lib/stb', 'lib/tinyexr', 'lib/glm', 'lib/glfw', 'lib/glad']
assert not git('diff', '--name-only', legacy_revision, source_revision, '--', *legacy_paths)
original_delta = git('diff', '--name-only', original_revision, source_revision, '--', *legacy_paths).decode().splitlines()
assert original_delta == ['src/Sirius.Core/Coordinate/PHCT002A.h']
starfield = 'src/Sirius.Render/Texture/Starfield.png'
assert git('rev-parse', original_revision + ':' + starfield).decode().strip() == entries[starfield]['git_blob'] == '2e652a5d1ba4ac1c64df640694d7584f30eb8a0f'
assert entries[starfield]['sha256'] == 'a87d68cb9e9e5ac2be1b7f1d5bdee2a8eabf07639d567643575102f202a1702c'
write(reference / 'legacy-origin.json', {'runtime_origin_revision': legacy_revision, 'build_scaffold_revision': source_revision,
      'byte_identical_path_domains': legacy_paths, 'original_revision': original_revision,
      'original_runtime_delta_paths': original_delta, 'delta': 'PHCT002A.h adds <algorithm>; not an untouched full original producer',
      'starfield': {'path': starfield, **entries[starfield]}})
for revision, path, destination in [
    ('70d156f', 'docs/SPECIFICATION.md', 'specification-70d156f.md'),
    ('8bd56ca', 'docs/SPECIFICATION.md', 'specification-8bd56ca.md')]:
    (reference / destination).write_bytes(git('show', revision + ':' + path))
(reference / 'named-e736ffe-transition.patch').write_bytes(git('diff', 'c24c0623f870106e4b1c5718e4b5b884fabae8db', 'e736ffe0c318dcd637f1453beb3523ad14a15047', '--', 'src/sirius/render'))
assert (reference / 'named-e736ffe-transition.patch').stat().st_size

dependency_heads = {'ftxui': 'cdf28903a7781f97ba94d30b79c3a4b0c97ccce7',
                    'json': '9cca280a4d0ccf0c08f47a99aa71d1b0e52f8d03',
                    'googletest': 'f8d7d77c06936315286eb55f8de22cd23c188571'}
dependencies = {}
dependency_files = []
for name, expected in dependency_heads.items():
    path = root / 'bin/linux-gcc/_deps' / (name + '-src')
    head = query(['git', '-C', str(path), 'rev-parse', 'HEAD']).decode().strip()
    status = query(['git', '-C', str(path), 'status', '--porcelain', '--untracked-files=normal']).decode()
    assert head == expected and not status, (name, head, status)
    records = {}
    for file in sorted(path.rglob('*')):
        relative = file.relative_to(path)
        if '.git' in relative.parts or file.is_dir():
            continue
        assert file.is_file() and not file.is_symlink(), file
        records[str(relative)] = identity(file)
        dependency_files.append(str(file))
    record = {'source_directory': str(path), 'git_head': head, 'git_status': status, 'files': records,
              'file_count': len(records), 'bytes': sum(e['bytes'] for e in records.values()), 'file_records_sha256': seal(records)}
    target = reference / 'dependency-manifests' / (name + '.json')
    write(target, record)
    dependencies[name] = {k: v for k, v in record.items() if k != 'files'}
    dependencies[name]['manifest'] = str(target)

defaults = {
    'render': {'width': 1920, 'height': 1080, 'samplesPerPixel': 64, 'tileSize': 64, 'threadCount': 0, 'outputPath': 'render.ppm'},
    'metric': {'name': 'Schwarzschild', 'mass': 1.0, 'spin': 0.0, 'charge': 0.0, 'lambda': 0.0,
               'temperatureModel': 'NovikovThorne', 'diskTemperature': 50000.0, 'throatRadius': 1.0,
               'warpVelocity': 0.5, 'bubbleRadius': 1.0, 'bubbleSigma': 0.5},
    'observer': {'distance': 50.0, 'inclination': 90.0, 'azimuth': 0.0, 'fov': 60.0},
    'postprocess': {'enableBloom': True, 'bloomIntensity': 0.3, 'bloomThreshold': 0.3, 'exposure': 1.0,
                    'contrast': 1.0, 'saturation': 1.0, 'tonemapper': 'ACES'},
    'backend': {'preferred': 'auto', 'enableDenoiser': False, 'cudaDevice': 0},
    'volumetric': {'enabled': False, 'hOverR': 0.1, 'hPower': 0.25, 'tauMidplane': 10.0, 'samples': 64,
                   'enableTurbulence': False, 'enableCorona': False},
    'film': {'enabled': False, 'preset': 'Interstellar', 'grainIntensity': 0.15, 'halationStrength': 0.15, 'vignetteStrength': 0.3}}
write(work / 'configs' / 'declared-defaults.json', defaults)
cases = []
for name, fmt, required in [('kerr-png', 'png', True), ('kerr-exr', 'exr', True), ('alcubierre-png', 'png', True), ('kerr-ppm-optional', 'ppm', False)]:
    kerr = name.startswith('kerr')
    width, height, requested, effective = (512, 512, 32, 25) if kerr else (1920, 1080, 64, 64)
    configs, render_argv = {}, {}
    for role, binary in [('legacy', 'Sirius'), ('ported', 'sirius')]:
        config = json.loads(json.dumps(defaults))
        config['render'].update(width=width, height=height, samplesPerPixel=requested,
                                outputPath=str(root / 'renders/forensic-reconstruction-b739f92' / role / (name + '.' + fmt)))
        config['metric'].update(name='Kerr' if kerr else 'Alcubierre', spin=0.9 if kerr else 0.0)
        config['observer'].update(distance=30.0 if kerr else 20.0, inclination=80.0, fov=60.0 if kerr else 90.0)
        config['backend']['preferred'] = 'cpu'
        path = work / 'configs' / role / (name + '.json')
        write(path, config)
        configs[role] = str(path)
        render_argv[role] = [str(build / 'renderers/bin' / binary), '--config', str(path), '--no-color', '--verbose', 'render', '--cpu']
    output = lambda role: json.loads(Path(configs[role]).read_text())['render']['outputPath']
    cases.append({'name': name, 'format': fmt, 'required': required, 'width': width, 'height': height,
                  'requested_spp': requested, 'effective_spp': effective, 'rays_per_renderer': width * height * effective,
                  'sample_rule': 'int(sqrt(float(requested_spp))) squared, deterministic subpixel grid',
                  'source_scene': original_revision + (':render_test.sh Kerr block' if kerr else ':render_demo.sh Alcubierre block'),
                  'format_scope': ('Kerr public PNG scene; EXR/PPM explicitly selected writer variants' if kerr else
                                   'Public Alcubierre script emits PPM; PNG is an explicit new-reconstruction writer variant'),
                  'configs': configs, 'future_render_argv_review_only': render_argv,
                  'future_compare_argv_after_owned_render_receipts': [str(build / 'comparator/image_compare'), 'compare', fmt, str(width), str(height), output('legacy'), output('ported')]})
assert sum(case['rays_per_renderer'] * 2 for case in cases if case['required']) == 291_635_200

preserved = root / 'attestations/scientific-core/4bcd987/p4-cpu-radiation-current-source/files/out/p4-cpu-radiation-current-source'
owner_hashes = {'run_owned.py': 'b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d',
                'owned_guard.py': '747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056'}
for name, expected in owner_hashes.items():
    assert identity(preserved / name)['sha256'] == expected
    shutil.copy2(preserved / name, work / name)
assert not git('status', '--porcelain', '--untracked-files=normal')
runner_revision = git('rev-parse', 'HEAD').decode().strip()
commands = {
    'configure-renderers': [['/usr/bin/cmake', '--preset', 'linux-gcc', '-S', str(source), '-B', str(build / 'renderers'),
                            '-DSIRIUS_BUILD_LEGACY=ON', '-DCMAKE_CUDA_COMPILER=NOTFOUND', '-DCMAKE_DISABLE_FIND_PACKAGE_Vulkan=TRUE',
                            '-DFETCHCONTENT_FULLY_DISCONNECTED=ON',
                            *['-DFETCHCONTENT_SOURCE_DIR_' + name.upper() + '=' + dependencies[name]['source_directory'] for name in dependencies]]],
    'build-renderers': [['/usr/bin/cmake', '--build', str(build / 'renderers'), '--target', 'Sirius', 'sirius', '--parallel', '2']],
    'build-comparator': [
        ['/usr/bin/gcc-14', '-O2', '-DNDEBUG', '-I' + str(source / 'lib/tinyexr'), '-c', str(source / 'lib/tinyexr/miniz.c'), '-o', str(build / 'comparator/miniz.o')],
        ['/usr/bin/g++-14', '-std=c++17', '-O2', '-DNDEBUG', '-fno-fast-math', '-ffp-contract=off', '-pthread',
         '-isystem', str(source / 'lib/stb'), '-isystem', str(source / 'lib/tinyexr'), str(work / 'image_compare.cpp'), str(build / 'comparator/miniz.o'), '-o', str(build / 'comparator/image_compare')]],
    'comparator-controls': [[str(build / 'comparator/image_compare'), 'fixtures', '{owned_folder}/fixtures'],
                            ['See phase.py: 17 comparisons, exact self/metadata/one-pixel outcomes and invalid-input refusal']]}
tools = {}
for name in ('/usr/bin/gcc-14', '/usr/bin/g++-14', '/usr/bin/cmake', '/usr/bin/ninja', '/usr/bin/git', '/usr/bin/python3.12'):
    path = Path(name).resolve(strict=True)
    tools[name] = {'resolved_path': str(path), **identity(path)}
manifest = {
    'evidence_kind': 'NEW_RECONSTRUCTION', 'prepared_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'status': 'PREPARED_NOT_EXECUTED', 'historical_tapes_recovered': False, 'historical_july_manifest_or_comparator_recovered': False,
    'claim_limit': 'New comparison of explicitly surviving legacy and early-ported CPU components only; no original tape identity, current-physics identity, scientific or release qualification.',
    'original_specification': '70d156f18fb8b72873df22cf0720fd922b0ad72c:docs/SPECIFICATION.md sections 1.6, 6(2), 7',
    'explicit_gate_replacement': '8bd56cadb4b1de5fd0f0f016dfec918062fd1b42; owner approval to relax original gate not recovered',
    'original_revision': original_revision, 'legacy_runtime_origin': legacy_revision, 'source_export_revision': source_revision,
    'runner_revision': runner_revision, 'runner_is_not_export_revision': True, 'source_manifest': str(work / 'source-manifest.json'),
    'legacy_origin_proof': str(reference / 'legacy-origin.json'), 'dependencies': dependencies, 'tools': tools,
    'codec_pins': {p: entries[p] for p in ['lib/stb/stb_image.h', 'lib/stb/stb_image_write.h', 'lib/tinyexr/tinyexr.h', 'lib/tinyexr/miniz.c', 'lib/tinyexr/miniz.h']},
    'starfield': {'path': starfield, **entries[starfield]}, 'runtime_cwd': str(source),
    'child_environment': {'set': {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'TZ': 'UTC', 'GIT_CEILING_DIRECTORIES': str(root / 'out'),
                                   'GIT_CONFIG_NOSYSTEM': '1', 'GIT_TERMINAL_PROMPT': '0', 'PYTHONDONTWRITEBYTECODE': '1'},
                          'unset': ['GIT_DIR', 'GIT_WORK_TREE', 'GIT_INDEX_FILE', 'GIT_OBJECT_DIRECTORY', 'GIT_ALTERNATE_OBJECT_DIRECTORIES'],
                          'reject': 'Nonempty SIRIUS_*, provider/GTest/CTest/loader overrides rejected by preserved owner; build flags/toolchain/include-path overrides rejected by phase.py. No HOME override.'},
    'git_discovery': 'Plain export has no .git; ceiling at canonical out prevents parent canonical Git identity from being used by source children. Source revision is explicit, never inferred from runner HEAD.',
    'configuration_layering': 'defaults -> explicit complete JSON -> SIRIUS environment -> CLI; reject SIRIUS overrides; CPU explicit in config and CLI',
    'worker_behavior': 'Both SessionConfig conversions leave threadCount at zero/auto rather than copying render.threadCount. Keep declared 0; observe and record actual worker count in any future owned render, not a promised fixed thread count.',
    'cases': cases, 'required_case_count': 3, 'required_render_count': 6, 'required_total_rays': 291_635_200,
    'optional_ppm_is_not_replacement': True, 'render_authorized_by_this_packet': False,
    'future_render_acceptance_dependencies': ['Root-reviewed separate one-shot ownership/deadlines before any of six required renders; no smaller substitute scene',
        'Actual built ELF/library/compile/cache provenance and binary hashes, complete config/environment/asset loaded without fallback',
        'Before/after input seals, worker observation, clean terminal ownership and successful complete writer output',
        'Record raw file SHA256 per role, then compare every declared pixel/channel using the validated comparator; preserve result and files',
        'Any differences require classification; no blanket expectation across e736ffe or later intentional physics changes'],
    'named_intentional_transition': {'before': 'c24c0623f870106e4b1c5718e4b5b884fabae8db', 'after': 'e736ffe0c318dcd637f1453beb3523ad14a15047',
        'behavior': 'Reduced-step geodesic rejection continues rather than terminating, attempts still count max_steps; Morris-Thorne Cartesian CPU embedding enabled',
        'not_in_selected_build_pair': True, 'no_new_measurement_claim': True},
    'comparator_contract': {'png': 'Full noninterlaced RGBA8, exact decoded 8-bit channels',
        'exr': 'Four A/B/G/R ZIP HALF scanline channels with declared origin/windows; full decoded RGBA IEEE754 binary32 words under sealed TinyEXR, no epsilon/finite substitution',
        'ppm': 'Full P6 RGB8, maxval255, exact dimensions/payload/EOF; optional legacy format, not a PNG/EXR substitute',
        'encoded_identity_separate': True, 'metadata_checked': 'Format/channel/compression/window/dimension domain plus bounded valid encoding; arbitrary ancillary metadata may differ and is not pixel identity',
        'writer_semantics': 'EXR stores linear RGBA HALF; PNG uses sRGB RGB/linear alpha rounding; PPM uses gamma2.2 RGB8. Compare within format; no claim about prequantization buffers.',
        'controls': 'Each format self and same-pixels/different-metadata accepted; valid one-represented-sample perturbation rejected with exactly one pixel/sample mismatch; missing/wrong dimensions refused; wrong format and truncated PNG refused',
        'limits': 'Ordinary nonsymlink input <=256MiB; dimensions<=16384 and <=64Mi samples, all requested scenes inside domain'},
    'preparation_commands': commands, 'preparation_limits': {'seconds_per_action': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152,
        'meaning': 'Imported operational diagnostic budget only; timeout is incomplete evidence, not a scientific mismatch or relaxation of the original Vulkan300s criterion'},
    'preserved_owner': {'source_directory': str(preserved), 'sha256': owner_hashes},
    'root_commands': [['/usr/bin/python3.12', str(work / 'run_action.py'), name] for name in commands],
    'no_agent_build_test_render_execution': True}
write(work / 'preparation.json', manifest)
packet_paths = [p for p in sorted(work.rglob('*')) if p.is_file()]
frozen_inputs = sorted(set([str(p) for p in packet_paths] + dependency_files + list(tools) +
                          [str(work / 'actions.json'), str(work / 'pins.json')]))
actions = {name: {'revision': runner_revision, 'kind': 'check', 'expected_exit': 0,
                 'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
                 'scope': 'NEW_RECONSTRUCTION ' + name + ' only; no original tapes, render, scientific or release qualification',
                 'argv': ['/usr/bin/python3.12', str(work / 'phase.py'), name, '{folder}'],
                 'frozen_inputs': frozen_inputs} for name in commands}
write(work / 'actions.json', actions)
pins = {str(p): identity(p) for p in sorted(work.rglob('*')) if p.is_file()}
pins.update({name: identity(Path(name)) for name in dependency_files + list(tools)})
write(work / 'pins.json', {'files': pins, 'file_count': len(pins), 'file_records_sha256': seal(pins)})
print(json.dumps({'prepared': True, 'export_files': len(entries), 'export_bytes': 43_614_536,
                  'preparation': identity(work / 'preparation.json'), 'pins': identity(work / 'pins.json'),
                  'dependencies': {name: dependencies[name]['file_records_sha256'] for name in dependencies}}, indent=2))
