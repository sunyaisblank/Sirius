"""Freeze the authenticated provider inputs for the original owned build guard."""
from pathlib import Path
import json, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent; S = W / 'upstream/mesa-25.2.8'
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
files = {str(p.relative_to(R)) for p in S.rglob('*') if p.is_file()}
inputs = files | {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
inputs |= {'/usr/bin/gcc-14', '/usr/bin/g++-14', '/usr/bin/ninja', '/usr/bin/cmake',
    '/usr/bin/ar', '/usr/bin/ld', '/usr/lib/llvm-21/bin/llvm-config',
    '/usr/lib/llvm-21/lib/libLLVM.so.21.1', str(R / 'bin/toolchains/python/bin/python3'),
    str(R / 'bin/toolchains/python/bin/meson')}
env = {'CC': '/usr/bin/gcc-14', 'CXX': '/usr/bin/g++-14',
    'PATH': str(R / 'bin/toolchains/python/bin') + ':/usr/lib/llvm-21/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
    'LLVM_CONFIG': '/usr/lib/llvm-21/bin/llvm-config', 'PYTHONDONTWRITEBYTECODE': '1'}
action = {'revision': revision, 'kind': 'build', 'expected_exit': 0,
    'limits': {'seconds': 1800, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'argv': ['/usr/bin/python3', '-B', str(W / 'build_provider.py')],
    'frozen_inputs': sorted(inputs), 'software_provider': False, 'build_environment': env,
    'scope': 'One private Mesa25.2.8/GCC14/LLVM21 provider build, exact existing IEEE32 prerequisite plus guarded exact acquire-release lowering. All original workgroup rendezvous, atomics and shader programs remain. No baseline/provider installation, source adoption, Sirius compilation or runtime/scientific/release acceptance.'}
(W / 'actions.json').write_text(json.dumps({'build': action}, indent=2) + '\n')
print(json.dumps({'revision': revision, 'immutable_upstream_files': len(files),
                  'new_provider_builds': 1, 'baseline_provider_builds': 0}))
