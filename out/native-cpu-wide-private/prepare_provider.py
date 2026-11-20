"""Apply the reviewed provider delta and recover its exact AX inputs."""

from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile

from resume_raw import BASELINE, R, check_space, pin, read_json, restore, write_json
from run_value_abi import expect_pin


HERE = Path(__file__).resolve().parent
QUALIFICATION = R / 'attestations/qualification'
MESA = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
AX = R / 'out/native-cpu-closed-all-provider/provider/source-contract.json'
ARCHIVES = (
    ('a47e219b/native-cpu-dispatch-admission-token-finite-compatibility/component.zip', range(3), ''),
    ('5ae04321/native-cpu-endpoint-header-frame-timeout/actual-products.zip', (3,), ''),
    ('86843614/native-cpu-local-identities-finite-compatibility/actual-products.zip', range(4, 7), 'files/'),
    ('363aa0c2/native-cpu-closed-all-original-c8/products.zip', (7,), ''),
)


def fingerprint_tree(tree):
    identities = {str(path.relative_to(tree)): pin(path) for path in tree.rglob('*') if path.is_file()}
    data = json.dumps(identities, sort_keys=True, separators=(',', ':')).encode()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def restore_member(bundle, member, output, record):
    expected = {key: record[key] for key in ('bytes', 'sha256')}
    if output.exists():
        expect_pin(output, expected)
    else:
        with bundle.open(member) as stream:
            restore(stream, output, expected)


def main():
    scratch, = (Path(arg).resolve() for arg in sys.argv[1:])
    if scratch == R / 'out' or not scratch.is_relative_to(R / 'out') or not scratch.is_dir():
        raise RuntimeError('Pass an existing per-task directory under out/')
    check_space()
    contract = read_json(AX)
    basis = fingerprint_tree(MESA)
    if basis != contract['full_candidate_source_identity']:
        raise RuntimeError('Mesa sources differ from the sealed AX tree')

    patch = HERE / 'provider.patch'
    names = [line.removeprefix('+++ b/') for line in patch.read_text().splitlines()
             if line.startswith('+++ b/')]
    candidate = scratch / 'mesa'
    shutil.copytree(MESA, candidate, symlinks=True, copy_function=os.link)
    for name in names:
        target = candidate / name
        detached = target.with_name(target.name + '.detached')
        shutil.copy2(target, detached)
        detached.replace(target)
    relative = str(candidate.relative_to(R))
    subprocess.run(['/usr/bin/git', 'apply', '--directory=' + relative, str(patch)], check=True)
    if fingerprint_tree(MESA) != basis:
        raise RuntimeError('The original Mesa sources changed')

    objects = contract['native_objects']
    for name, indexes, prefix in ARCHIVES:
        with zipfile.ZipFile(QUALIFICATION / name) as bundle:
            for index in indexes:
                record = objects[index]
                path = record['path']
                restore_member(bundle, prefix + path, R / path, record)

    header_inputs = read_json(BASELINE)['replay']['source_inputs']
    header_prefix = 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch/'
    archive = QUALIFICATION / '536ac050/native-cpu-bridge-source-and-boundaries/actual-native-products.zip'
    include = scratch / 'programs'
    with zipfile.ZipFile(archive) as bundle:
        for stage in ('transport', 'endpoint', 'dopri_phase'):
            name = stage + '.program.hpp'
            path = header_prefix + name
            restore_member(bundle, path, include / name, header_inputs[path])

    adapter = scratch / 'admitted_rows.o'
    compile_argv = ['/usr/bin/g++-14', '-std=c++2c', '-O2', '-fPIC', '-march=x86-64', '-mtune=znver4',
                    '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-isystem', str(include), '-c',
                    str(HERE / 'admitted_rows.cpp'), '-o', str(adapter)]
    subprocess.run(compile_argv, check=True, env=dict(os.environ, TMPDIR=str(scratch)))
    write_json(scratch / 'provider-inputs.json', {'basis': basis, 'candidate': fingerprint_tree(candidate),
        'provider_patch': pin(patch), 'original_objects': objects, 'adapter': pin(adapter),
        'adapter_source': pin(HERE / 'admitted_rows.cpp'), 'compile_argv': compile_argv,
        'changed_files': {name: pin(candidate / name) for name in names}})


if __name__ == '__main__':
    main()
