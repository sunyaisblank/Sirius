"""Restore the verified SIMD16 stage objects without restoring compiler toolchains."""

from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

from replay_native import EMIT, SEALED
from resume_raw import R, check_space, pin, read_json, restore, write_json
from run_value_abi import MATH_CAPSULE, expect_pin
from prepare_provider import AX

OBJECTS = {
    'adapter-object': 'wide-math.o',
    'emit-transport': 'transport.o',
    'emit-endpoint': 'endpoint.o',
    'emit-dopri_phase': 'dopri_phase.o',
}


def require_repo_directory(path, root, description):
    if not path.is_dir() or path == root or not path.is_relative_to(root):
        raise RuntimeError(f'Use an existing {description} directory below {root}')


def materialize(archive, member, path, expected):
    if path.exists():
        expect_pin(path, expected)
        return
    temporary = path.with_suffix('.restore')
    if temporary.exists():
        raise RuntimeError(f'Pending restore exists: {temporary}')
    try:
        with archive.open(member) as stream:
            restore(stream, temporary, {key: expected[key] for key in ('bytes', 'sha256')})
        temporary.rename(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def restore_objects(capsule, scratch):
    check_space()
    require_repo_directory(capsule, R / 'attestations', 'qualification capsule')
    require_repo_directory(scratch, R / 'out', 'task scratch')
    evidence = read_json(capsule / 'MANIFEST.json')
    archived = read_json(SEALED)
    expect_pin(SEALED, evidence['sealed_proof'])
    archive_path = capsule / 'provider-objects.zip'
    expect_pin(archive_path, evidence['archive'])
    expect_pin(Path('/usr/bin/ld'), evidence['linker'])

    inputs = evidence['native_replay']['results']
    results = {}
    outputs = {}
    with ZipFile(archive_path) as archive:
        if archive.namelist() != list(OBJECTS.values()):
            raise RuntimeError('Unexpected native capsule members or order')
        for action, name in OBJECTS.items():
            frozen = archived['results'][action]['output']
            if any(inputs[action]['output'][key] != frozen[key] for key in ('bytes', 'sha256')):
                raise RuntimeError(f'Invalid sealed pin for {action}')
            output = (R / frozen['path']) if action == 'adapter-object' else scratch / 'native' / name
            materialize(archive, name, output, frozen)
            results[action] = {**inputs[action],
                               'output': {'path': str(output.relative_to(R)), **pin(output)}}
            outputs[action] = output

    core = read_json(AX)['closed_dp_math_object']
    with ZipFile(MATH_CAPSULE / 'products.zip') as archive:
        materialize(archive, core['path'], R / core['path'], core)

    linked = scratch / 'native/simd16-component.o'
    argv = ['/usr/bin/ld', '-r', str(outputs['adapter-object']), str(R / core['path']),
            *(str(outputs[action]) for action in EMIT)]
    pending = linked.with_suffix('.pending.o')
    if pending.exists():
        raise RuntimeError(f'Pending component link exists: {pending}')
    frozen_link = archived['results']['link']['output']
    try:
        subprocess.run([*argv, '-o', str(pending)], check=True)
        expect_pin(pending, frozen_link)
        if linked.exists():
            expect_pin(linked, frozen_link)
        else:
            pending.rename(linked)
    finally:
        if pending.exists():
            pending.unlink()
    results['link'] = {**inputs['link'],
                       'output': {'path': str(linked.relative_to(R)), **pin(linked)}}
    receipt = scratch / 'native-replay.json'
    write_json(receipt, {**evidence['native_replay'], 'results': results,
                        'restored_from': {'path': str(archive_path.relative_to(R)), **pin(archive_path)},
                        'restore_linker_argv': [*argv, '-o', str(pending)]})
    print(f'Verified and restored {len(OBJECTS)} sealed inputs into {receipt.relative_to(R)}')


if __name__ == '__main__':
    capsule, scratch = (Path(arg).resolve() for arg in sys.argv[1:])
    restore_objects(capsule, scratch)
