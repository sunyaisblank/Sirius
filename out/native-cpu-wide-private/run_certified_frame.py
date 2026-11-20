"""Run the original cold frame with the sealed native-owned renderer and certified continuation."""

from pathlib import Path
import json
import re
import shutil
import subprocess
import sys

from run_original_c8 import (GUARD, LIMITS, ORIGINAL_READER, QUALIFICATION, R, check_space,
                             expect_pin, pin, read_json, write_json)
from run_original_controls import AX_FRAME
from run_original_frame import prepare_resources


SEALED_SCENE = QUALIFICATION / '2ec9eea2/native-cpu-simd16-original-gates/MANIFEST.json'


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
    reference_scene = read_json(SEALED_SCENE)['original_cold_frame']['scene']
    candidate_path = scratch / 'certified-frame-binary.json'
    candidate = read_json(candidate_path)
    executable = scratch / 'frame-assets/sirius_render_tests'
    identity = candidate['execution_image']
    assert candidate['revision'] == revision and R / identity['path'] == executable
    expect_pin(executable, identity)
    for file in candidate['original_workload']:
        expect_pin(R / file['path'], file)
    historical['frame_executable'] = identity
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
    frozen = [executable, provider, icd, profile, ORIGINAL_READER, AX_FRAME, SEALED_SCENE, Path(__file__).resolve(),
              scratch / 'all8-acceptance.json', candidate_path, *(resources / name for name in files)]
    argv = [original['argv'][0], '-DSIRIUS_RENDER_TESTS=' + str(executable), *original['argv'][2:]]
    write_json(work / 'actions.json', {'frame': {'revision': revision, 'kind': 'test', 'expected_exit': 0,
        'limits': LIMITS, 'argv': argv, 'software_provider': True, 'disable_software_shader_cache': True,
        'runtime_icd': str(icd), 'runtime_provider': str(provider), 'frozen_inputs': list(map(str, frozen)),
        'build_environment': {'TMPDIR': str(temporary_render), 'SIRIUS_CERTIFIED_SPHERICAL_ESCAPE': '1'},
        'scope': 'Original 128x128 cold Kerr frame and unchanged limits with certified spherical continuation'}})
    subprocess.run([sys.executable, '-B', str(work / 'run_owned.py'), 'frame'], check=False)

    folder = work / 'frame'
    owner = read_json(folder / 'owner.json')
    stdout = folder / 'stdout.log'
    stderr = folder / 'stderr.log'
    xml = folder / 'gtest.xml'
    text = stdout.read_text()
    error_text = stderr.read_text()
    diagnostics = error_text.splitlines()
    publication = [line for line in diagnostics if line.startswith('[Vulkan] Certified spherical escapes:')]
    problems = []
    render_evidence, scenes = [], []
    try:
        if (owner['status'] != 'terminal' or not owner['source_unchanged'] or
                not owner['whole_inputs_unchanged'] or not owner['provider_mappings'] or
                owner['remaining_owned_processes'] or owner['cleanup_errors']):
            raise RuntimeError('The original owner did not verify inputs, process closure and provider mapping')
        render_evidence = [json.loads(line.split('Render evidence: ', 1)[1])
                           for line in text.splitlines() if line.startswith('[Vulkan] Render evidence: ')]
        scenes = [json.loads(line.split('Scene evidence: ', 1)[1]) for line in text.splitlines()
                  if line.startswith('[Session] Scene evidence: ')]
        if (owner['passed'] and not scenes) or any(scene != reference_scene for scene in scenes):
            raise RuntimeError('The original Kerr frame scene was not verified')
        if owner['passed']:
            frame, = render_evidence
            if (frame['metric'], frame['width'], frame['height'], frame['precision'], frame['route'],
                    frame['work_items']) != ('Kerr', 128, 128, 'fp32', 'retained', 16384):
                raise RuntimeError('The original retained route or tile count was not verified')
            if frame['accepted_intervals'] <= 0 or not any(
                    re.fullmatch(r'\[Vulkan\] Certified spherical escapes: [1-9]\d*', line)
                    for line in publication):
                raise RuntimeError('The adaptive and certified continuation routes were not both observed')
            sys.path.insert(0, str(ORIGINAL_READER.parent))
            import native_evidence as reader
            reader.CASES['frame'] = 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'
            reader.xml_properties(xml, 'frame')
    except (AssertionError, KeyError, OSError, RuntimeError, ValueError) as exc:
        problems.append(f'{type(exc).__name__}: {exc}')
    accepted = owner['passed'] and not problems
    progress = [int(value) for value in re.findall(r'\bTile (\d+)/16384\b', text)]
    checkpoints = [i for i, line in enumerate(diagnostics) if line.startswith('[Vulkan] Retained progress: ')]
    last_snapshot = diagnostics[checkpoints[-1]:] if checkpoints else []
    snapshot_prefixes = ('[Vulkan] Retained progress:', '[Vulkan] Retained cumulative physical stage:',
                         '[Vulkan] Retained cumulative shared Endpoint/Dense', '[Vulkan] Retained stage prefixes:')
    result = {'source': revision, 'accepted': accepted, 'evidence_errors': problems, 'provider': all8['provider'],
              'executable': pin(executable), 'resources': files, 'original_all8': pin(scratch / 'all8-acceptance.json'),
              'owner': pin(folder / 'owner.json'), 'stdout': pin(stdout), 'stderr': pin(stderr),
              'xml': pin(xml) if xml.exists() else None,
              'stop_reason': owner['stop_reason'], 'child_exit': owner['child_exit'], 'wall_seconds': owner['wall_seconds'],
              'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'], 'limits': LIMITS,
              'last_reported_tiles': 16384 if accepted else max(progress, default=0),
              'total_tiles': 16384, 'scene': scenes[0] if scenes else None,
              'render_evidence': render_evidence,
              'candidate_preflight': pin(candidate_path),
              'certified_publication': publication,
              'last_stage_snapshot': [line for line in last_snapshot if line.startswith(snapshot_prefixes)]}
    write_json(scratch / 'cold-frame-result.json', result)
    shutil.rmtree(temporary_render)
    return result


if __name__ == '__main__':
    location, = (Path(arg).resolve() for arg in sys.argv[1:])
    if location == R / 'out' or not location.is_relative_to(R / 'out'):
        raise RuntimeError('Use the working directory for this provider build under out/')
    sys.exit(0 if run(location)['accepted'] else 1)
