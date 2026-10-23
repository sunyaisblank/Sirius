from pathlib import Path
import hashlib
import json
import shlex
import shutil
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
output = root / 'bin/linux-gcc/lavapipe-exact-pure-outlining/outline_tool'
output.parent.mkdir(parents=True, exist_ok=True)
compiler = Path(shutil.which('g++-14')).resolve(strict=True)
config = Path('/usr/lib/llvm-21/bin/llvm-config')
assert subprocess.check_output([str(config), '--version'], text=True).strip() == '21.1.8'
flags = shlex.split(subprocess.check_output([str(config), '--cxxflags'], text=True))
selected = []
for flag in flags:
    if flag.startswith('-I'):
        selected.extend(['-isystem', flag[2:]])
    else:
        selected.append(flag)
links = shlex.split(subprocess.check_output(
    [str(config), '--ldflags', '--libs', '--system-libs'], text=True))
argv = [str(compiler), *selected, '-std=c++2c', '-O2', '-Wall', '-Wextra',
        '-Wpedantic', '-Werror', str(work / 'outline_tool.cpp'),
        str(work / 'verify_inverse.cpp'),
        '-o', str(output), *links]
pin = lambda p: {'bytes': p.stat().st_size,
                 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
record = {'argv': argv, 'compiler': pin(compiler), 'llvm_config': pin(config),
          'scope': 'One private offline LLVM transformer executable; no Sirius, shader, observer, provider or numerical test compilation. Primary source/tool pins, not an exhaustive operating-system toolchain qualification.'}
result = subprocess.run(argv)
record['exit'] = result.returncode
if result.returncode == 0:
    record['product'] = {'path': str(output), **pin(output)}
(work / 'build-invocation.json').write_text(json.dumps(record, indent=2) + '\n')
raise SystemExit(result.returncode)
