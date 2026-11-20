"""Build the SIMD16 provider under the original Mesa build limits."""

from pathlib import Path
import os
import shutil
import subprocess
import sys

from prepare_provider import AX, R, fingerprint_tree
from resume_raw import GUARD, check_space, pin, read_json, write_json
from run_value_abi import expect_pin


PREFIX = R / 'bin/linux-gcc'
BUILD = PREFIX / 'lavapipe-acqrel-ordering'
TARGET = 'src/gallium/targets/lavapipe/libvulkan_lvp.so'
LIMITS = {'seconds': 1800, 'rss_kib': 4194304, 'reserve_kib': 2097152}


def build(scratch):
    check_space()
    if subprocess.check_output(['git', 'status', '--porcelain']).strip():
        raise RuntimeError('Run from the clean candidate checkout')
    prepared = read_json(scratch / 'provider-inputs.json')
    component = read_json(scratch / 'provider-component.json')
    candidate = scratch / 'mesa'
    assert prepared['candidate'] == fingerprint_tree(candidate)
    files = [R / item['path'] for item in read_json(AX)['native_objects']]
    for path, item in zip(files, prepared['original_objects'], strict=True):
        expect_pin(path, item)
    wide = R / component['output']['path']
    expect_pin(wide, component['output'])
    env = {'CC': '/usr/bin/gcc-14', 'CXX': '/usr/bin/g++-14', 'LLVM_CONFIG': '/usr/lib/llvm-21/bin/llvm-config',
           'PATH': '/usr/lib/llvm-21/bin:' + os.environ['PATH'], 'PYTHONDONTWRITEBYTECODE': '1',
           'PYTHONPATH': str(scratch / 'meson-tool'), 'TMPDIR': str(scratch)}
    meson = ['/usr/bin/python3', '-B', '-m', 'mesonbuild.mesonmain']
    assert subprocess.check_output([*meson, '--version'], env=dict(os.environ, **env), text=True).strip() == '1.9.1'
    config = [*meson, 'setup', str(BUILD), str(candidate),
              '-Dvulkan-drivers=swrast', '-Dgallium-drivers=llvmpipe', '-Dplatforms=', '-Dglx=disabled',
              '-Degl=disabled', '-Dgbm=disabled', '-Dopengl=false', '-Dgles1=disabled', '-Dgles2=disabled',
              '-Dllvm=enabled', '-Dshared-llvm=enabled', '-Dllvm-orcjit=false', '-Dbuild-tests=false',
              '-Dmicrosoft-clc=disabled', '-Dbuildtype=release', '-Db_ndebug=true',
              '-Dc_args=-DLP_VULKAN_IEEE32=1', '-Dprefix=' + str(BUILD / 'uninstalled-prefix')]
    work = scratch / 'provider-build'
    work.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', work / 'run_owned.py')
    (work / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    inputs = [*files, wide, scratch / 'provider-inputs.json', scratch / 'provider-component.json',
              *(path for path in candidate.rglob('*') if path.is_file()),
              *(path for path in (scratch / 'meson-tool').rglob('*') if path.is_file()),
              *(Path(name) for name in (env['CC'], env['CXX'], env['LLVM_CONFIG'], '/usr/bin/ninja', '/usr/bin/ld'))]
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    common = {'revision': revision, 'kind': 'build', 'expected_exit': 0, 'software_provider': False,
              'limits': LIMITS, 'build_environment': env, 'frozen_inputs': sorted(map(str, inputs)),
              'scope': 'Private SIMD16 Vulkan provider build using the original Mesa limits and configuration'}
    actions = {}
    if not (BUILD / 'build.ninja').exists():
        actions['configure'] = dict(common, argv=config)
    else:
        info = read_json(BUILD / 'meson-info/meson-info.json')
        assert Path(info['directories']['source']).resolve() == candidate
    actions['build'] = dict(common, argv=['/usr/bin/ninja', '-C', str(BUILD), '-j2', TARGET])
    write_json(work / 'actions.json', actions)
    for action in actions:
        subprocess.run([sys.executable, '-B', str(work / 'run_owned.py'), action], check=True)
    write_json(scratch / 'provider-build.json', {'source': revision, 'provider': pin(BUILD / TARGET),
        'build_owner': pin(work / 'build/owner.json'), 'candidate_source': prepared['candidate'],
        'configuration': config, 'build_limits': LIMITS})


if __name__ == '__main__':
    location, = (Path(arg).resolve() for arg in sys.argv[1:])
    if location == R / 'out' or not location.is_relative_to(R / 'out'):
        raise RuntimeError('Use the working directory for this provider build under out/')
    build(location)
