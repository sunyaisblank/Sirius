# Derived from the existing cache/configuration planner; no new native compiler action.
"""Seal unchanged build products and configuration/tool/source inputs before Meson."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
B = R / 'bin/linux-gcc/lavapipe-acqrel-ordering'
P = W / 'provider'
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
review = json.loads((P / 'review-admission.json').read_text())
assert review['passed']
base = json.loads((R / 'out/native-cpu-masked-arithmetic/provider/actions.json').read_text())['build']
inputs = {name for name in base['frozen_inputs'] if not (R / name).is_relative_to(B)}
inputs |= set(json.loads((P / 'source-contract.json').read_text())['source_inputs'])
inputs |= {str(path.relative_to(R)) for path in S.rglob('*') if path.is_file()}
inputs |= {str(path.relative_to(R)) for path in W.rglob('*') if path.is_file()
           and 'build' not in path.relative_to(W).parts}
mutable_directories = {'meson-private', 'meson-logs', 'meson-info', 'meson-uninstalled'}
mutable_names = {'build.ninja', 'compile_commands.json', '.ninja_log', '.ninja_deps', '.ninja_lock'}
reused = []
for path in B.rglob('*'):
    if not path.is_file():
        continue
    relative = path.relative_to(B)
    if relative.parts[0] in mutable_directories or path.name in mutable_names:
        continue
    inputs.add(str(path.relative_to(R)))
    reused.append(str(path.relative_to(R)))
for item in json.loads((P / 'source-contract.json').read_text())['native_objects']:
    inputs.add(item['path'])
python_root = R / 'bin/toolchains/python'
inputs |= {str(path.relative_to(R)) for path in
           (python_root / 'lib/python3.12/site-packages/mesonbuild').rglob('*')
           if path.is_file() and '__pycache__' not in path.parts}
inputs |= {str(python_root / 'bin/python3'), str(python_root / 'bin/meson'),
           '/usr/bin/ninja', '/usr/bin/pkg-config', '/usr/bin/cmake',
           '/usr/lib/llvm-21/bin/llvm-config', '/usr/lib/llvm-21/lib/libLLVM.so.21.1'}
assert all((R / name).is_file() for name in inputs)
env = {'CC': '/usr/bin/gcc-14', 'CXX': '/usr/bin/g++-14',
       'PATH': str(python_root / 'bin') + ':/usr/lib/llvm-21/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
       'LLVM_CONFIG': '/usr/lib/llvm-21/bin/llvm-config', 'PYTHONDONTWRITEBYTECODE': '1'}
action = {
    'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'kind': 'configure_only', 'expected_exit': 0,
    'limits': {'seconds': 60, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'argv': ['/usr/bin/python3', '-B', str(W / 'configure_provider.py')],
    'frozen_inputs': sorted(inputs), 'software_provider': False, 'build_environment': env,
    'scope': 'Reconfigure existing private Mesa build once to reuse8 native objects with three raw SIMD objects calling the checked dense C++ arithmetic plus the combined original adapter/helper object; same existing build options/toolchain. All prior object/archive/provider contents remain sealed. No shader/native execution or benefit/science/release qualification.',
}
assert not (P / 'actions.json').exists()
(P / 'actions.json').write_text(json.dumps({'configure': action}, indent=2) + '\n')
for name in ('run_owned.py', 'owned_guard.py'):
    (P / name).write_text((W / name).read_text())
print(json.dumps({'sealed_inputs': len(inputs), 'reused_build_files': len(reused),
                  'new_configuration_options': 0}))
