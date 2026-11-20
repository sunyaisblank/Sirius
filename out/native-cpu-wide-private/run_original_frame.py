"""Run the original cold 128x128 Kerr frame once after the original all8."""

from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile

from run_original_c8 import (BINARY_ARCHIVE, GUARD, LIMITS, ORIGINAL_READER, QUALIFICATION,
                             R, check_space, expect_pin, pin, read_json, restore, sealed, write_json)
from run_original_controls import AX_FRAME, AX_ROOT


RESOURCE_ARCHIVE = QUALIFICATION / '34ef39fb/native-cpu-arithmetic-sharing-finite-outcome/new-actual-products-01.zip'
RESOURCE_ROOT = 'bin/linux-gcc/native-cpu-arithmetic-sharing/frame/resources/'
MODULE_ARCHIVE = QUALIFICATION / 'bfa095ba/native-cpu-row-ownership-finite-outcome/new-actual-products-01.zip'


def prepare_resources(scratch, record):
    frame = scratch / 'frame-assets'
    frame.mkdir(exist_ok=True)
    executable = frame / 'sirius_render_tests'
    exe = record['frame_executable']
    if executable.exists():
        expect_pin(executable, exe)
    else:
        with zipfile.ZipFile(BINARY_ARCHIVE) as archive, archive.open('files/' + exe['path']) as stream:
            restore(stream, executable, {key: exe[key] for key in ('bytes', 'sha256')})
        executable.chmod(0o755)
    module_records = record['final_frame_result']['module_joins']
    linked = executable.read_bytes()
    with zipfile.ZipFile(MODULE_ARCHIVE) as archive:
        for program in module_records.values():
            data = archive.read('sha256/' + program['sha256'])
            assert len(data) == program['bytes'] and hashlib.sha256(data).hexdigest() == program['sha256']
            assert linked.count(data) == 1

    bundle = record['bundle']
    expect_pin(R / bundle['path'], bundle)
    with zipfile.ZipFile(R / bundle['path']) as archive:
        files = json.loads(sealed(archive, bundle['members'], AX_ROOT + '/cold-frame/resource-restore.json'))['files']
        original_actions = json.loads(sealed(archive, bundle['members'], AX_ROOT + '/cold-frame/actions.json'))
    resources = frame / 'resources'
    seen = {}
    with zipfile.ZipFile(RESOURCE_ARCHIVE) as archive:
        for name, expected in files.items():
            target = resources / name
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                expect_pin(target, expected)
                seen[expected['sha256']] = target
                continue
            if name == 'assets/Starfield.png':
                asset = R / name
                expect_pin(asset, expected)
                os.link(asset, target)
            elif expected['sha256'] in seen:
                os.link(seen[expected['sha256']], target)
            else:
                with archive.open(RESOURCE_ROOT + name) as stream:
                    restore(stream, target, expected)
            expect_pin(target, expected)
            seen[expected['sha256']] = target
    return executable, resources, files, original_actions['frame']


def run(scratch):
    check_space()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run the original frame from the clean C8 candidate checkout')
    all8 = read_json(scratch / 'all8-acceptance.json')
    c8 = read_json(scratch / 'c8-acceptance.json')
    expect_pin(scratch / 'c8-acceptance.json', all8['original_c8'])
    expect_pin(scratch / 'abba/benefit.json', all8['original_abba'])
    expect_pin(scratch / 'abba/decimal-agreement.json', all8['decimal'])
    assert c8['source'] == all8['source'] == revision and all8['all8_original_controls_pass']
    assert len(all8['remaining7']) == 7
    provider = R / all8['provider']['path']
    expect_pin(provider, all8['provider'])
    icd = scratch / 'lvp_icd.json'
    historical = read_json(AX_FRAME)
    executable, resources, files, original = prepare_resources(scratch, historical)
    profile = R / 'tests/operational/vulkan_runtime_profile.cmake'
    assert original['argv'][0] == '/usr/bin/cmake' and original['argv'][2:] == ['-P', str(profile)]
    assert original['limits'] == LIMITS and original['expected_exit'] == 0
    temporary_render = R / 'renders' / scratch.name
    temporary_render.mkdir()
    work = scratch / 'cold-frame'
    work.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', work / 'run_owned.py')
    (work / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    frozen = [executable, provider, icd, profile, ORIGINAL_READER, AX_FRAME, Path(__file__).resolve(),
              scratch / 'all8-acceptance.json', *(resources / name for name in files)]
    argv = [original['argv'][0], '-DSIRIUS_RENDER_TESTS=' + str(executable), *original['argv'][2:]]
    write_json(work / 'actions.json', {'frame': {'revision': revision, 'kind': 'test', 'expected_exit': 0,
        'limits': LIMITS, 'argv': argv, 'software_provider': True, 'disable_software_shader_cache': True,
        'runtime_icd': str(icd), 'runtime_provider': str(provider), 'frozen_inputs': list(map(str, frozen)),
        'build_environment': {'TMPDIR': str(temporary_render)},
        'scope': 'Original cold software Vulkan frame, 128x128 Kerr one SPP with the original assertions and limits'}})
    subprocess.run([sys.executable, '-B', str(work / 'run_owned.py'), 'frame'], check=False)

    folder = work / 'frame'
    owner = read_json(folder / 'owner.json')
    stdout = folder / 'stdout.log'
    stderr = folder / 'stderr.log'
    xml = folder / 'gtest.xml'
    assert owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['provider_mappings'] and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    text = stdout.read_text()
    scenes = [json.loads(line.split('Scene evidence: ', 1)[1]) for line in text.splitlines()
              if line.startswith('[Session] Scene evidence: ')]
    assert scenes and all(scene['metric'] == 'Kerr' and scene['backend'] == 'Vulkan' and
                         (scene['width'], scene['height'], scene['samples_per_pixel']) == (128, 128, 1)
                         for scene in scenes)
    if owner['passed']:
        sys.path.insert(0, str(ORIGINAL_READER.parent))
        import native_evidence as reader
        reader.CASES['frame'] = 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'
        reader.xml_properties(xml, 'frame')
    progress = [int(value) for value in re.findall(r'\bTile (\d+)/16384\b', text)]
    diagnostics = stderr.read_text().splitlines()
    checkpoints = [i for i, line in enumerate(diagnostics) if line.startswith('[Vulkan] Retained progress: ')]
    last_snapshot = diagnostics[checkpoints[-1]:] if checkpoints else []
    snapshot_prefixes = ('[Vulkan] Retained progress:', '[Vulkan] Retained cumulative physical stage:',
                         '[Vulkan] Retained cumulative shared Endpoint/Dense', '[Vulkan] Retained stage prefixes:')
    result = {'source': revision, 'accepted': owner['passed'], 'provider': all8['provider'],
              'executable': pin(executable), 'resources': files, 'original_all8': pin(scratch / 'all8-acceptance.json'),
              'owner': pin(folder / 'owner.json'), 'stdout': pin(stdout), 'stderr': pin(stderr),
              'xml': pin(xml) if xml.exists() else None,
              'stop_reason': owner['stop_reason'], 'child_exit': owner['child_exit'], 'wall_seconds': owner['wall_seconds'],
              'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'], 'limits': LIMITS,
              'last_reported_tiles': max(progress, default=0), 'total_tiles': 16384, 'scene': scenes[0],
              'last_stage_snapshot': [line for line in last_snapshot if line.startswith(snapshot_prefixes)]}
    write_json(scratch / 'cold-frame-result.json', result)
    shutil.rmtree(temporary_render)
    return result


if __name__ == '__main__':
    location, = (Path(arg).resolve() for arg in sys.argv[1:])
    if location == R / 'out' or not location.is_relative_to(R / 'out'):
        raise RuntimeError('Use the working directory for this provider build under out/')
    sys.exit(0 if run(location)['accepted'] else 1)
