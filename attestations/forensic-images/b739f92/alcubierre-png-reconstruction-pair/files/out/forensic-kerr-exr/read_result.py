"""Read one owned full-size reconstruction Kerr result; --self-decode is root-only CPU validation."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parent.parent
CASE = 'legacy-kerr-exr-original-size'
FOLDER = WORK / CASE


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load(path):
    return json.loads(path.read_text())


def identity(path):
    require(path.is_file(), 'missing ordinary file: ' + str(path))
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            size += len(block)
            digest.update(block)
    return {'bytes': size, 'sha256': digest.hexdigest()}


def path_identity(path):
    return {'resolved_path': str(path.resolve(strict=True)), **identity(path)}


def proc_birth(pid):
    try:
        stat = Path(f'/proc/{pid}/stat').read_text()
        return stat[stat.rfind(')') + 2:].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None


def map_record(line):
    match = re.fullmatch(r'([0-9a-f]+)-([0-9a-f]+)\s+([-rwxps]{4})\s+([0-9a-f]+)\s+([0-9a-f]+):([0-9a-f]+)\s+(\d+)\s+(.+)', line)
    require(match is not None, 'unparsed file mapping: ' + line)
    return {'permissions': match[3], 'device': os.makedev(int(match[5], 16), int(match[6], 16)),
            'inode': int(match[7]), 'path': match[8]}


def main():
    global CASE, FOLDER
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', help='owned case directory; defaults to the sole action in actions.json')
    parser.add_argument('--self-decode', action='store_true', help='root explicitly permits one bounded validated-comparator CPU self-decode')
    args = parser.parse_args()
    preparation = load(WORK / 'preparation.json')
    actions = load(WORK / 'actions.json')
    require(args.case is not None or len(actions) == 1, 'choose one declared --case')
    CASE = args.case if args.case is not None else next(iter(actions))
    require(CASE in actions and Path(CASE).name == CASE and CASE not in ('.', '..'), 'case must name one declared owned directory')
    FOLDER = WORK / CASE
    role = actions[CASE]['role']
    require(role in ('legacy', 'ported'), 'unsupported reconstruction producer role')
    action = actions[CASE]
    owner = load(FOLDER / 'owner.json')
    result = {'evidence_kind': 'NEW_RECONSTRUCTION', 'case': CASE, 'role': role,
              'classification': 'UNVERIFIED', 'complete_original_size_render': False,
              'historical_tape_identity_claimed': False, 'pair_pixel_identity_claimed': False,
              'scientific_or_release_qualification_claimed': False, 'proof_gaps': [],
              'reader_sha256': identity(Path(__file__).resolve())['sha256']}
    require(owner['status'] == 'terminal', 'owner is not terminal; do not inspect or decode a live output')
    require(owner['case'] == CASE and owner['limits'] == preparation['limits'] == action['limits'] ==
            {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}, 'wrong owner/case/limits')
    require(owner['source'] == {'revision': preparation['runner_revision'], 'status': ''}, 'wrong runner source identity')
    require(owner['source_unchanged'] and owner['whole_inputs_unchanged'], 'owner input/source change')
    require(not owner['cleanup_errors'] and owner['remaining_owned_processes'] == [], 'ownership retirement incomplete')
    births = owner['observed_births']
    require(births and all(item['absent'] and item['current_start_ticks'] != item['saved_start_ticks'] for item in births), 'observed birth not retired')
    require(Path('/proc/sys/kernel/random/boot_id').read_text().strip() == owner['boot_id'], 'fresh birth check requires the recorded boot')
    for item in births:
        require(proc_birth(item['pid']) != item['saved_start_ticks'], 'recorded process birth is still present')
    before, after = load(FOLDER / 'inputs-before.json'), load(FOLDER / 'inputs-after.json')
    require(before == after and before, 'whole input seals differ')
    current = {}
    for name, saved in before.items():
        path = Path(name) if Path(name).is_absolute() else ROOT / name
        resolved = str(path.resolve(strict=True))
        if resolved not in current:
            current[resolved] = identity(path)
        require(saved == {'resolved_path': resolved, **current[resolved]}, 'sealed input now differs: ' + name)
    pins = load(WORK / 'pins.json')['files']
    for name, saved in pins.items():
        path = Path(name)
        key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else name
        require(key in before and {k: before[key][k] for k in ('bytes', 'sha256')} == saved, 'pin absent from execution seals: ' + name)
    for name in action['frozen_inputs']:
        path = Path(name)
        key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else name
        require(key in before, 'declared frozen input not sealed: ' + name)
    expected_argv = [arg.replace('{folder}', str(FOLDER)) for arg in action['argv']]
    require(owner['argv'] == expected_argv, 'executed owner argv differs')
    result['ownership'] = {'terminal': True, 'clean_retirement': True, 'observed_births_absent': len(births),
                           'absence_scope': 'Recorded process births only; no global consumer-absence claim',
                           'input_seals_verified': len(before), 'pins_verified': len(pins),
                           'stop_reason': owner['stop_reason'], 'child_exit': owner['child_exit'],
                           'wall_seconds': owner['wall_seconds'], 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'],
                           'minimum_sampled_available_kib': owner['minimum_sampled_available_kib']}
    samples = [json.loads(line) for line in (FOLDER / 'samples.jsonl').read_text().splitlines()]
    require(samples and owner['peak_sampled_rss_kib'] >= max(item['rss_kib'] for item in samples) and
            owner['minimum_sampled_available_kib'] <= min(item['available_kib'] for item in samples), 'resource samples disagree with owner summary')
    if owner['stop_reason'] is None:
        require(owner['peak_sampled_rss_kib'] <= action['limits']['rss_kib'] and
                owner['minimum_sampled_available_kib'] >= action['limits']['reserve_kib'], 'accepted owner resource bound differs')
    result['ownership']['raw_sample_count'] = len(samples)
    config_path = Path(preparation['case']['configs'][role])
    config = load(config_path)
    require(preparation['case']['name'] == 'kerr-exr' and preparation['case']['format'] == 'exr', 'wrong reconstruction domain')
    require(config['render']['width'] == config['render']['height'] == 512 and config['render']['samplesPerPixel'] == 32 and
            config['render']['tileSize'] == 64 and config['render']['threadCount'] == 0 and
            config['backend']['preferred'] == 'cpu', 'original render domain differs')
    require(config['metric']['name'] == 'Kerr' and config['metric']['mass'] == 1 and config['metric']['spin'] == 0.9 and
            config['observer'] == {'distance': 30.0, 'inclination': 80.0, 'azimuth': 0.0, 'fov': 60.0}, 'scene configuration differs')
    packet = Path(preparation['source_packet'])
    source_manifest = load(packet / 'source-manifest.json')
    starfield = Path(preparation['cwd']) / 'src/Sirius.Render/Texture/Starfield.png'
    require(identity(starfield) == {k: source_manifest['files']['src/Sirius.Render/Texture/Starfield.png'][k]
                                  for k in ('bytes', 'sha256')}, 'wrong historical asset bytes')
    producer = preparation['producers'][role]
    elf = Path(producer['path'])
    require(identity(elf) == {k: producer[k] for k in ('bytes', 'sha256')}, 'wrong renderer bytes')
    bindings = load(FOLDER / 'live-bindings.json')
    require(len(bindings) == 1, 'expected one captured renderer binding')
    binding = bindings[0]
    for field in ('pid', 'start_ticks', 'pgrp', 'session'):
        require(binding['before'][field] == binding['after'][field], 'unstable live binding: ' + field)
    pid = binding['before']['pid']
    require(binding['before']['session'] == owner['child_pid'], 'renderer outside owned session')
    require(any(b['pid'] == pid and b['saved_start_ticks'] == binding['before']['start_ticks'] for b in births), 'renderer birth not joined to owner')
    require(binding['exe'] == str(elf) and binding['exe_bytes'] == producer['bytes'] and binding['exe_sha256'] == producer['sha256'], 'live ELF binding differs')
    require(binding['cwd'] == preparation['cwd'] and binding['argv'] == preparation['case']['future_render_argv_review_only'][role], 'live cwd/argv differs')
    raw_maps = (FOLDER / 'live-maps.txt').read_text().splitlines()
    elf_maps = [map_record(line) for line in raw_maps if line.split()[-1:] == [str(elf)]]
    require(elf_maps and any('x' in item['permissions'] for item in elf_maps), 'actual executable mapping absent')
    for item in elf_maps:
        require(item['inode'] == elf.stat().st_ino and item['device'] == elf.stat().st_dev, 'ELF map backing identity differs')
    ldd_path = Path(producer['ldd'])
    ldd_key = str(ldd_path.relative_to(ROOT)) if ldd_path.is_relative_to(ROOT) else str(ldd_path)
    require(ldd_key in before, 'role-specific LDD record was not frozen for execution')
    role_libraries = set()
    for line in ldd_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if re.fullmatch(r'linux-vdso\.so\.\d+ \(0x[0-9a-f]+\)', line):
            continue  # Kernel-provided mapping, not a backing-file library.
        match = re.fullmatch(r'(?:\S+\s+=>\s+)?(/[^\s()]+)\s+\(0x[0-9a-f]+\)', line)
        require(match is not None, 'unresolved or unparsed role LDD record: ' + line)
        role_libraries.add(str(Path(match[1]).resolve(strict=True)))
    require(role_libraries and role_libraries <= set(preparation['dynamic_libraries']), 'role dependency absent from frozen union')
    require(set(binding['libraries']) == role_libraries, 'mapped libraries differ from the producer-specific LDD set')
    executable_files = set()
    for line in raw_maps:
        fields = line.split(maxsplit=5)
        if len(fields) == 6 and fields[5].startswith('/') and 'x' in fields[1]:
            executable_files.add(map_record(line)['path'])
    require(executable_files == role_libraries | {str(elf)}, 'unexpected or missing executable file mapping')
    for name in sorted(role_libraries):
        saved = preparation['dynamic_libraries'][name]
        path, live = Path(name), binding['libraries'][name]
        require(identity(path) == saved == live['file'], 'mapped library bytes differ: ' + name)
        require(live['stat_inode'] == path.stat().st_ino and live['stat_dev'] == path.stat().st_dev, 'mapped library stat differs: ' + name)
        require(live['mappings'] and all(line in raw_maps for line in live['mappings']), 'library mapping record absent: ' + name)
        parsed = [map_record(line) for line in live['mappings']]
        require(any('x' in item['permissions'] for item in parsed), 'executable library mapping absent: ' + name)
        require(all(item['path'] == name and item['inode'] == live['stat_inode'] and item['device'] == live['stat_dev'] for item in parsed), 'library map backing differs: ' + name)
    trace = (FOLDER / 'runtime.trace').read_text()
    process_lines = [line for line in trace.splitlines() if re.match(rf'^{pid}\s', line)]
    exec_lines = [line for line in process_lines if 'execve(' + json.dumps(str(elf)) + ',' in line and re.search(r'\) = 0$', line)]
    require(len(exec_lines) == 1, 'successful actual renderer exec not proved')

    def successful_open(path):
        return [line for line in process_lines if 'openat(' in line and ', ' + json.dumps(path) + ',' in line and re.search(r'\) = \d+(?:<[^>]*>)?$', line)]

    config_opens = successful_open(str(config_path))
    asset_opens = successful_open('src/Sirius.Render/Texture/Starfield.png')
    require(config_opens and all('O_RDONLY' in line for line in config_opens), 'sealed config open not proved')
    require(asset_opens and all('O_RDONLY' in line and 'AT_FDCWD<' + preparation['cwd'] + '>' in line for line in asset_opens), 'sealed source-relative asset open not proved')
    out = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', (FOLDER / 'stdout.log').read_text())
    err = (FOLDER / 'stderr.log').read_text()
    require('[Session] Loaded starfield texture: 6000x3000 RGBA' in out, 'successful full asset decode not logged')
    require('Background will render as faint grey' not in out + err and 'Could not load starfield texture from any path' not in out + err, 'actual grey fallback occurred')
    workers = re.findall(r'Multi-threaded rendering enabled: (\d+) threads', out)
    started = re.findall(r'Starting (\d+) worker threads', out)
    require(len(workers) == len(started) == 1 and workers == started and int(workers[0]) > 0, 'worker count not joined')
    require('--cpu' in binding['argv'], 'forced CPU argument absent')
    if role == 'legacy':
        require('[Session] GPU disabled by config, using CPU rendering' in out, 'legacy forced CPU path not logged')
    require(re.search(r'Resolution:\s+512 x 512', out) and re.search(r'Samples:\s+32(?: spp)?', out) and
            re.search(r'Tiles:\s+64 \(spiral order\)', out) and re.search(r'Metric:\s+Kerr', out), 'runtime scene summary differs')
    tiles = [(int(n), int(total)) for n, total in re.findall(r'Tile (\d+)/(\d+)', out)]
    require(tiles and all(total == 64 and 0 < n <= 64 for n, total in tiles), 'invalid or absent tile progress')
    require([n for n, _ in tiles] == sorted(set(n for n, _ in tiles)), 'tile completion log is not strictly increasing')
    result['runtime_join'] = {'renderer_pid': pid, 'renderer_birth': binding['before']['start_ticks'], 'producer': producer,
        'config': {'path': str(config_path), **identity(config_path)}, 'asset': {'path': str(starfield), **identity(starfield)},
        'mapped_libraries': {name: preparation['dynamic_libraries'][name] for name in sorted(role_libraries)},
        'role_ldd_record': {'path': str(ldd_path), **identity(ldd_path)},
        'other_role_frozen_libraries_not_required_to_be_mapped': sorted(set(preparation['dynamic_libraries']) - role_libraries),
        'workers': int(workers[0]),
        'requested_spp': 32, 'effective_spp_from_sealed_source': 25, 'tiles_reported': [n for n, _ in tiles],
        'completed_tile_ray_count': max(n for n, _ in tiles) * 4096 * 25,
        'source_scope': 'Declared original CPU semantics; temperature options remain unwired and selected azimuth0 remains fixed',
        'binding_scope': 'Stable-birth live ELF/backing-library snapshot, sealed inputs, successful config/asset opens and source-derived CPU/sampler join; modern source has no separate CPU-selection log',
        'expected_failed_asset_lookups_allowed': ['Sirius.Render/Texture/Starfield.png', '../Sirius.Render/Texture/Starfield.png']}
    require(identity(WORK / 'run_owned.py') == identity(packet / 'run_owned.py'), 'owner differs from the preserved environment-rejection source')
    require(expected_argv[:3] == ['/usr/bin/env', '-C', preparation['cwd']], 'environment wrapper/cwd differs')
    required_environment = {'LC_ALL': 'C', 'TZ': 'UTC', 'PATH': '/usr/bin:/bin',
                            'GIT_CEILING_DIRECTORIES': str(ROOT / 'out'), 'GIT_CONFIG_NOSYSTEM': '1'}
    assignments = expected_argv[3:expected_argv.index('/usr/bin/strace')]
    require(all('=' in value for value in assignments), 'unparsed environment-wrapper argument')
    declared_environment = dict(value.split('=', 1) for value in assignments)
    require(declared_environment == required_environment and owner['declared_output_environment'] == {'GTEST_OUTPUT': ''}, 'declared relevant environment differs')
    observed_environment = binding.get('observed_relevant_environment')
    if observed_environment is not None:
        require(all(observed_environment.get(name) == value for name, value in required_environment.items()), 'live relevant environment assignment differs')
        prohibited = [name for name, value in observed_environment.items() if value and
                      (name.startswith(('SIRIUS_', 'VK_', 'LP_', 'MESA_', 'GALLIVM_', 'GTEST_', 'CTEST_')) or
                       name in {'LD_PRELOAD', 'LD_AUDIT', 'LD_LIBRARY_PATH'})]
        require(not prohibited, 'unexpected live config/provider/loader/test override: ' + ', '.join(sorted(prohibited)))
    result['runtime_join']['environment'] = {
        'relevant_assignments': required_environment,
        'config_provider_loader_test_overrides': 'Nonempty relevant overrides rejected by sealed preserved parent before constructing the child environment; GTEST_OUTPUT then set empty',
        'live_scoped_snapshot_present': observed_environment is not None,
        'scope': 'Renderer/config/provider/loader/test and declared locale/path/Git selectors; unrelated variables and credentials are outside this claim'}
    output = Path(action['output'])
    require(str(output) == config['render']['outputPath'], 'output path differs from opened config')
    if owner['stop_reason'] is not None:
        require(not owner['passed'], 'guard-stopped owner marked passed')
        classification = 'INCOMPLETE_DIAGNOSTIC_LIMIT' if owner['stop_reason'].startswith('diagnostic_') else 'INCOMPLETE_INTERRUPTED'
        result.update(classification=classification,
                      scientific_assertion_verdict='unavailable', output=({'path': str(output), **identity(output)} if output.exists() else None))
        result['proof_gaps'].append('No accepted terminal original-size output; guard stop is not a scientific mismatch or pixel comparison verdict')
        return result, 3
    require(owner['passed'] and owner['child_exit'] == 0, 'render child failed without an accepted terminal owner')
    require([n for n, _ in tiles] == list(range(1, 65)) and '[Session] All worker threads completed' in out, 'full original tile completion not proved')
    require('[Session] Wrote: ' + str(output) in out and 'Render complete: ' + str(output) in out, 'successful writer/terminal completion not logged')
    output_opens = successful_open(str(output))
    require(output_opens and any(all(flag in line for flag in ('O_WRONLY', 'O_CREAT', 'O_TRUNC')) for line in output_opens), 'actual output write open not proved')
    require(output.is_file() and not output.is_symlink(), 'ordinary EXR output absent')
    saved_output = identity(output)
    result['output'] = {'path': str(output), **saved_output}
    comparator = Path(preparation['comparator']['path'])
    build_receipt = load(packet / 'build-comparator/receipt.json')
    controls = load(packet / 'comparator-controls/receipt.json')
    controls_owner = load(packet / 'comparator-controls/owner.json')
    expected_binary = next(item for item in build_receipt['products'] if item['path'] == preparation['comparator']['copied_from'])
    require(identity(comparator) == {k: expected_binary[k] for k in ('bytes', 'sha256')}, 'validated comparator bytes differ')
    require(controls['passed'] and controls_owner['passed'] and controls_owner['status'] == 'terminal' and controls_owner['child_exit'] == 0 and
            controls_owner['source_unchanged'] and controls_owner['whole_inputs_unchanged'] and not controls_owner['cleanup_errors'] and
            controls_owner['remaining_owned_processes'] == [] and controls_owner['stop_reason'] is None, 'comparator controls not accepted')
    argv = [str(comparator), 'compare', 'exr', '512', '512', str(output), str(output)]
    result['self_decode'] = {'argv': argv, 'binary': identity(comparator), 'passed': False,
                            'meaning': 'Full-file decode/domain validation only; self-comparison is not legacy/ported or historical tape identity'}
    if not args.self_decode:
        result['classification'] = 'COMPLETE_RENDER_DECODE_PENDING'
        result['proof_gaps'].append('Root must run this reader with --self-decode to validate the full ZIP/HALF RGBA EXR decoded to exact binary32 words')
        return result, 2
    start = time.monotonic()
    decoded = subprocess.run(argv, cwd=ROOT, env={'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'TZ': 'UTC'}, capture_output=True, text=True, timeout=300)
    result['self_decode'].update(exit=decoded.returncode, stdout=decoded.stdout, stderr=decoded.stderr,
                                 wall_seconds=time.monotonic() - start)
    require(identity(output) == saved_output and identity(comparator) == {k: expected_binary[k] for k in ('bytes', 'sha256')}, 'decode input/tool changed')
    require(decoded.returncode == 0 and not decoded.stderr, 'EXR self-decode refused')
    detail = json.loads(decoded.stdout)
    require(detail == {'equal': True, 'format': 'exr', 'sample_domain': 'decoded-RGBA-IEEE754-binary32-words-from-HALF', 'width': 512, 'height': 512,
                       'channels': 4, 'compared_pixels': 262144, 'compared_samples': 1048576,
                       'mismatched_pixels': 0, 'mismatched_samples': 0, 'encoded_files_equal': True}, 'unexpected full EXR decode domain/result')
    result['self_decode']['passed'] = True
    result.update(classification='COMPLETE_ORIGINAL_SIZE_SINGLE_RENDER', complete_original_size_render=True,
                  interpretation='One declared-role reconstruction EXR completed and decoded as exact binary32 words without epsilon; this reader establishes no result for another render or any cross-producer comparison')
    return result, 0


if __name__ == '__main__':
    try:
        payload, code = main()
    except Exception as error:
        payload, code = {'evidence_kind': 'NEW_RECONSTRUCTION', 'case': CASE, 'classification': 'EVIDENCE_GAP',
                         'complete_original_size_render': False, 'historical_tape_identity_claimed': False,
                         'pair_pixel_identity_claimed': False, 'scientific_or_release_qualification_claimed': False,
                         'proof_gaps': [str(error)]}, 1
    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(code)
