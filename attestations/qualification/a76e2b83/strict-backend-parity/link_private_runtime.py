"""Link the strict Linux render test with the sealed native-owned modules."""

import json
from pathlib import Path
import shlex
import subprocess
import sys
from zipfile import ZipFile


ROOT = Path.cwd()
WORK = ROOT / 'out/strict-backend-science'
PRESET = ROOT / 'bin/linux-ci'
SELECTOR = WORK / 'native-provider/selector'
OUTPUT = WORK / 'native-provider/matched-runtime'
TARGET = 'src/sirius/app/sirius_render_tests'
PRIVATE_STAGES = {'transport', 'endpoint', 'dopri_phase'}
ORIGINAL_FILES = ('tests/render/app_render_execution_test.cpp',
                  'tests/operational/vulkan_runtime_profile.cmake')
CASE = 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'

sys.path.insert(0, str(ROOT / 'out/native-cpu-wide-private'))
import prepare_owned_frame_modules as selector
from resume_raw import pin, read_json, write_json
from run_original_frame import MODULE_ARCHIVE
from run_value_abi import expect_pin


def link_command():
    query = ['/usr/bin/ninja', '-C', str(PRESET), '-t', 'commands', TARGET]
    matches = []
    for line in subprocess.check_output(query, text=True).splitlines():
        args = shlex.split(line)
        if '-o' in args and args[args.index('-o') + 1] == TARGET:
            matches.append(args)
    command, = matches
    if command[:2] != [':', '&&']:
        raise RuntimeError('Unexpected strict renderer link recipe')
    return command[2:command.index('&&', 2)]


def link():
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain']):
        raise RuntimeError('Expected an exact clean strict candidate')
    for path in ORIGINAL_FILES:
        original = subprocess.check_output(['git', 'show', f'2ec9eea2:{path}'])
        if original != (ROOT / path).read_bytes():
            raise RuntimeError(f'The original workload changed: {path}')

    selector.PRESET = PRESET
    selector.HEADER = PRESET / 'src/sirius/backend/retained/retained_kernels.h'
    SELECTOR.mkdir(exist_ok=True)
    proof = SELECTOR / 'owned-frame-modules.json'
    if not proof.exists():
        selector.prepare(SELECTOR)
    selection = read_json(proof)
    if selection['revision'] != revision:
        raise RuntimeError('The selector was compiled from a different source')
    compute = SELECTOR / 'retained_compute.o'
    expect_pin(compute, selection['object'])

    policy = read_json(WORK / 'certified-link/link.json')
    renderer = ROOT / policy['object']['path']
    expect_pin(renderer, policy['object'])
    if policy['checkout'] != revision:
        raise RuntimeError('The continuation policy was compiled from a different source')
    OUTPUT.mkdir(exist_ok=True)
    execution = OUTPUT / 'sirius_render_tests'
    command = link_command()
    output_index = command.index('-o')
    command[output_index + 1] = str(execution)
    command[output_index:output_index] = [str(compute), str(renderer)]
    subprocess.run(command, cwd=PRESET, check=True)

    executable = execution.read_bytes()
    if b'SIRIUS_CERTIFIED_SPHERICAL_ESCAPE\x00' in executable:
        raise RuntimeError('Expected the strict build without the development environment switch')
    with ZipFile(MODULE_ARCHIVE) as archive:
        for stage, record in selection['original_renderer_modules'].items():
            data = (archive.read('sha256/' + record['sha256']) if stage in PRIVATE_STAGES else
                    (selector.HEADER.parent / f'retained_{stage}.spv').read_bytes())
            if not selector.matches(data, record) or executable.count(data) != 1:
                raise RuntimeError(f'Could not bind the original private selection for {stage}')
    listing = subprocess.check_output([execution, '--gtest_filter=' + CASE, '--gtest_list_tests'],
                                      text=True)
    if not any(line.strip() == CASE.split('.', 1)[1] for line in listing.splitlines()):
        raise RuntimeError('The original cold Vulkan case is missing')

    resources = OUTPUT / 'resources'
    if not resources.exists():
        subprocess.run(['/usr/bin/cp', '-al', str((PRESET / TARGET).parent / 'resources'),
                        str(resources)], check=True)
    if not (resources / 'assets/Starfield.png').samefile(ROOT / 'assets/Starfield.png'):
        raise RuntimeError('Expected the existing Starfield asset inode')
    write_json(OUTPUT / 'link.json', {
        'revision': revision, 'scope': 'Private provider runtime diagnostic with strict render objects',
        'executable': {'path': str(execution.relative_to(ROOT)), **pin(execution)},
        'selector': selection, 'render_policy': policy['object'],
        'link_command': command, 'original_case': CASE,
        'original_workload_inputs': {path: pin(ROOT / path) for path in ORIGINAL_FILES},
        'register_output': listing,
    })


if __name__ == '__main__':
    link()
