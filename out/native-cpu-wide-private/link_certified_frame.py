"""Relink the original render test case with the current code and sealed native modules."""

from pathlib import Path
import shlex
import subprocess
import sys
from zipfile import ZipFile

from prepare_owned_frame_modules import HEADER, PRESET, SIMD16, matches
from run_original_frame import MODULE_ARCHIVE
from resume_raw import R, check_space, pin, read_json, write_json
from run_value_abi import expect_pin

TARGET = 'src/sirius/app/sirius_render_tests'
CASE = 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'
TEST_SOURCE = 'tests/render/app_render_execution_test.cpp'
PROFILE = 'tests/operational/vulkan_runtime_profile.cmake'
SEALED_FRAME_REVISION = '2ec9eea2'


def identity(path):
    return {'path': str(path.relative_to(R)), **pin(path)}


def compile_candidate(scratch, owned):
    check_space()
    if any(n == R / 'out' or not n.is_relative_to(R / 'out') or not n.is_dir()
           for n in (scratch, owned)):
        raise RuntimeError('Use existing per-task directories below out/')
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Relink from a clean checkout')
    for name in (TEST_SOURCE, PROFILE):
        frozen = subprocess.check_output(['git', 'show', f'{SEALED_FRAME_REVISION}:{name}'])
        if frozen != (R / name).read_bytes():
            raise RuntimeError(f'Original frame workload has changed: {name}')

    module_proof = owned / 'owned-frame-modules.json'
    module_receipt = read_json(module_proof)
    if module_receipt['revision'] != revision:
        raise RuntimeError('Relinked shader selector does not share the checkout revision')
    selector = owned / 'retained_compute.o'
    expect_pin(selector, module_receipt['object'])

    query = ['/usr/bin/ninja', '-C', str(PRESET), '-t', 'commands', TARGET]
    choices = []
    for line in subprocess.check_output(query, text=True).splitlines():
        args = shlex.split(line)
        if '-o' in args and args[args.index('-o') + 1] == TARGET:
            choices.append(args)
    link, = choices
    if link[:2] != [':', '&&'] or not Path(link[2]).is_absolute():
        raise RuntimeError('Unexpected renderer test link recipe')
    argv = link[2:link.index('&&', 2)]
    if any(part in {';', '|', '||', ':', '>'} for part in argv):
        raise RuntimeError('Cannot replay a compound link recipe')
    assets = scratch / 'frame-assets'
    assets.mkdir(exist_ok=True)
    executable = assets / 'sirius_render_tests'
    if executable.exists():
        raise RuntimeError('Frame executable already exists')
    output = argv.index('-o')
    argv[output + 1] = str(executable)
    argv.insert(output, str(selector))
    subprocess.run(argv, cwd=PRESET, check=True)

    copy = executable.read_bytes()
    if b'the required runtime profile cannot admit retained Kerr rendering' not in copy:
        raise RuntimeError('The original required Vulkan assertion is missing')
    expected = read_json(SIMD16)['original_c8']['selected_modules']
    owned_stages = {'transport', 'endpoint', 'dopri_phase'}
    with ZipFile(MODULE_ARCHIVE) as archive:
        for stage, record in expected.items():
            data = (archive.read('sha256/' + record['sha256']) if stage in owned_stages
                    else (HEADER.parent / f'retained_{stage}.spv').read_bytes())
            if not matches(data, record) or copy.count(data) != 1:
                raise RuntimeError(f'Original renderer module join failed: {stage}')
    if b'SIRIUS_CERTIFIED_SPHERICAL_ESCAPE\0' not in copy:
        raise RuntimeError('The development continuation selector is missing')
    argv_list = [str(executable), '--gtest_filter=' + CASE, '--gtest_list_tests']
    listed = subprocess.check_output(argv_list, text=True)
    if not any(line.strip() == CASE.split('.', 1)[1] for line in listed.splitlines()):
        raise RuntimeError('The original single render case is unavailable')
    write_json(scratch / 'certified-frame-binary.json', {
        'revision': revision, 'execution_image': identity(executable),
        'original_workload': [identity(R / name) for name in (TEST_SOURCE, PROFILE)],
        'link_query': query, 'link_argv': argv, 'link_directory': str(PRESET),
        'selector': identity(selector), 'selector_receipt': identity(module_proof),
        'embedded_modules': expected, 'list_argv': argv_list, 'list_output': listed,
    })


if __name__ == '__main__':
    scratch, selector = (Path(arg).resolve() for arg in sys.argv[1:])
    compile_candidate(scratch, selector)
