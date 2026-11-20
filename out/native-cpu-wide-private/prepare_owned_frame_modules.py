"""Prepare the exact archived native-owned module selector for the cold-frame candidate."""

from pathlib import Path
import hashlib
import json
import shlex
import struct
import subprocess
import sys
from zipfile import ZipFile

from resume_raw import R, check_space, pin, read_json, write_json
from run_original_c8 import ORIGINAL
from run_original_frame import AX_FRAME, MODULE_ARCHIVE

PRESET = R / 'bin/linux-gcc'
HEADER = PRESET / 'src/sirius/backend/retained/retained_kernels.h'
SELECTOR = R / 'out/native-cpu-row-ownership/host-selector.patch'
SOURCE = R / 'src/sirius/backend/retained_compute.cpp'
SIMD16 = R / 'attestations/qualification/2ec9eea2/native-cpu-simd16-original-gates/MANIFEST.json'


def matches(data, record):
    return len(data) == int(record['bytes']) and hashlib.sha256(data).hexdigest() == record['sha256']


def prepare(scratch):
    check_space()
    if scratch == R / 'out' or not scratch.is_relative_to(R / 'out') or not scratch.is_dir():
        raise RuntimeError('Use an existing per-task directory below out/')
    historical = read_json(AX_FRAME)['final_frame_result']['module_joins']
    selected = read_json(SIMD16)['original_c8']['selected_modules']
    inputs = {}
    named = {'Transport': 'transport', 'Endpoint': 'endpoint', 'DopriPhase': 'dopri_phase'}
    other_stages = ('camera', 'dense', 'initialize', 'ray_camera')
    for stage in other_stages:
        path = HEADER.parent / f'retained_{stage}.spv'
        data = path.read_bytes()
        if not matches(data, selected[stage]):
            raise RuntimeError(f'The ordinary {stage} module differs from the pinned frame input')
        inputs[stage] = data
    with ZipFile(MODULE_ARCHIVE) as archive:
        for display, stage in named.items():
            row = historical[display]
            data = archive.read('sha256/' + row['sha256'])
            if not matches(data, row) or not matches(data, selected[stage]):
                raise RuntimeError(f'The archived native-owned {stage} module pin differs')
            inputs[stage] = data

    generated = scratch / 'retained_kernels.h'
    with generated.open('x') as output:
        output.write(f'#pragma once\n#include "{HEADER}"\nnamespace sirius::backend::retained_program {{\n')
        for display, stage in named.items():
            data = inputs[stage]
            words = struct.unpack('<' + 'I' * (len(data) // 4), data)
            output.write(f'inline constexpr std::array<std::uint32_t, {len(words)}> k{display}NativeOwnedShader{{{{\n')
            for i in range(0, len(words), 16):
                output.write(','.join(str(value) + 'u' for value in words[i:i+16]) + ',\n')
            output.write('}};\n')
        output.write('}\n')

    patched = scratch / SOURCE.name
    if patched.exists():
        raise RuntimeError(f'An owned selector already exists: {patched}')
    patch = ['/usr/bin/patch', '--fuzz=0', '--batch', '--output=' + str(patched),
             str(SOURCE), str(SELECTOR)]
    subprocess.run(patch, check=True)
    roster = json.loads((PRESET / 'compile_commands.json').read_text())
    recipe, = [row for row in roster if row['file'] == str(SOURCE)]
    argv = shlex.split(recipe['command'])
    argv[argv.index(str(SOURCE))] = str(patched)
    argv.insert(1, '-I' + str(scratch))
    artifact = scratch / 'retained_compute.o'
    argv[argv.index('-o') + 1] = str(artifact)
    subprocess.run(argv, cwd=recipe['directory'], check=True)
    complete = artifact.read_bytes()
    for name, data in inputs.items():
        if complete.count(data) != 1:
            raise RuntimeError(f'Expected one complete embedded copy of {name}')
    write_json(scratch / 'owned-frame-modules.json', {
        'scope': 'Compile the archived native-owned module selector with the current build inputs',
        'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'selector': pin(SELECTOR), 'source': pin(SOURCE), 'generated_header': pin(generated),
        'patch_argv': patch, 'compile_argv': argv, 'object': pin(artifact),
        'c8_authority': pin(ORIGINAL), 'module_authority': pin(AX_FRAME),
        'original_renderer_modules': {name: selected[name] for name in sorted(inputs)},
    })


if __name__ == '__main__':
    scratch, = (Path(arg).resolve() for arg in sys.argv[1:])
    prepare(scratch)
