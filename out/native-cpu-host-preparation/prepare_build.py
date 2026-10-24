"""Seal one C++ compile, one-member archive replacement and two original links."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
cache = root / 'bin/linux-gcc'
products = cache / 'native-cpu-host-preparation'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (work / 'prepared-build.json').exists()
assert not (work / 'host-build-result.json').exists()
assert not (work / 'build/actions.json').exists()
source = doc(work / 'host-source-proof.json')
assert source['full_inverse_restores_matched_source'] and source['all_other_source_bytes_unchanged']
assert pin(root / source['derived_path']) == source['derived']
authority = doc(work / 'reuse-authority.json')
inputs = dict(source['source_inputs'])
inputs[source['derived_path']] = source['derived']
for name, expected in authority['source_inputs'].items():
    assert pin(root / name) == expected, name
    inputs[name] = expected
original = doc(root / authority['matched_receipt'])
frame = doc(root / authority['frame_recipe'])
assert original['passed'] and original['new_cpp_compilations'] == 1 and original['observer_compilations'] == 0
commands = [dict(item) for item in original['commands']]
assert len(commands) == 2
object_path = products / 'matched/retained_compute.cpp.o'
archive_path = products / 'matched/libsirius_backend.a'
observer_path = products / 'reused/audited-observer.o'
host_path = products / 'matched/audited_backend_tests'
frame_path = products / 'frame/sirius_render_tests'
compile_args = list(commands[0]['argv'])
compile_args[compile_args.index('-c') + 1] = str(root / source['derived_path'])
compile_args[compile_args.index('-o') + 1] = str(object_path)
old_archive = str(cache / 'native-cpu-row-ownership/matched/libsirius_backend.a')
old_observer = str(cache / 'native-cpu-row-ownership/matched/audited-observer.o')
host_args = list(commands[1]['argv'])
assert host_args.count(old_archive) == 2 and host_args.count(old_observer) == 1
host_args = [str(archive_path) if arg == old_archive else
             str(observer_path) if arg == old_observer else arg for arg in host_args]
host_args[host_args.index('-o') + 1] = str(host_path)
frame_args = list(frame['argv'])
assert frame_args.count(old_archive) == 1
frame_args = [str(archive_path) if arg == old_archive else arg for arg in frame_args]
frame_args[frame_args.index('-o') + 1] = str(frame_path)
inverse = list(compile_args)
inverse[inverse.index('-c') + 1] = commands[0]['argv'][commands[0]['argv'].index('-c') + 1]
inverse[inverse.index('-o') + 1] = commands[0]['argv'][commands[0]['argv'].index('-o') + 1]
assert inverse == commands[0]['argv']
inverse = [old_archive if arg == str(archive_path) else
           old_observer if arg == str(observer_path) else arg for arg in host_args]
inverse[inverse.index('-o') + 1] = commands[1]['argv'][commands[1]['argv'].index('-o') + 1]
assert inverse == commands[1]['argv']
inverse = [old_archive if arg == str(archive_path) else arg for arg in frame_args]
inverse[inverse.index('-o') + 1] = frame['argv'][frame['argv'].index('-o') + 1]
assert inverse == frame['argv']
planned = [{'argv': compile_args, 'cwd': commands[0]['cwd']},
           {'argv': ['/usr/bin/ar', 'rDs', str(archive_path), str(object_path)], 'cwd': str(cache)},
           {'argv': host_args, 'cwd': commands[1]['cwd']},
           {'argv': frame_args, 'cwd': frame['cwd']}]
for path in (object_path, archive_path, observer_path, host_path, frame_path):
    assert not path.exists(), path
for name in ('matched', 'reused', 'frame'):
    (products / name).mkdir(parents=True, exist_ok=True)
for path in work.rglob('*'):
    if path.is_file():
        inputs[str(path.relative_to(root))] = pin(path)
# Record actual cached dependency names without invoking a compiler or build.
dependencies = subprocess.check_output(['/usr/bin/ninja', '-C', str(cache), '-t', 'deps',
    'src/sirius/backend/CMakeFiles/sirius_backend.dir/retained_compute.cpp.o'], text=True)
header_names = [line.strip() for line in dependencies.splitlines()[1:] if line.strip()]
assert len(header_names) >= 300
header_names += [str(path) for path in (root / 'out/native-cpu-row-ownership/headers').iterdir() if path.is_file()]
for name in header_names:
    path = Path(name)
    if not path.is_absolute():
        path = cache / path
    assert path.is_file(), path
    inputs[str(path.relative_to(root)) if path.is_relative_to(root) else str(path)] = pin(path)
# Freeze the exact compiler/binutils implementation, start files and libraries,
# as well as its dynamically resolved dependencies. No flag/tool selection.
tools = [Path(name) for name in ('/usr/bin/g++-14', '/usr/bin/ar', '/usr/bin/ranlib',
    '/usr/bin/as', '/usr/bin/ld', '/usr/bin/ninja', '/usr/bin/python3',
    '/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus')]
for path in (Path('/usr/lib/gcc/x86_64-linux-gnu/14')).rglob('*'):
    if path.is_file():
        tools.append(path)
for path in tools:
    resolved = path.resolve(strict=True)
    inputs[str(resolved)] = pin(resolved)
    if resolved.stat().st_mode & 0o111:
        result = subprocess.run(['/usr/bin/ldd', str(resolved)], capture_output=True, text=True)
        for name in re.findall(r'(/[^\s()]+)', result.stdout):
            dependency = Path(name)
            if dependency.is_file():
                resolved_dependency = dependency.resolve(strict=True)
                inputs[str(resolved_dependency)] = pin(resolved_dependency)
assert inputs[str(Path('/usr/bin/g++-14').resolve())]['sha256'] == '916ef3c03c2d024857af2e084d3648204033b9c79e0e62575974f567f99adc84'
assert inputs[str(Path('/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus').resolve())]['sha256'] == '52fda0dfd1df95d8d9370c39de67e08b6e37dcfe25a2ef977f8451521159e834'
inputs['/etc/ld.so.cache'] = pin(Path('/etc/ld.so.cache'))
for name, expected in inputs.items():
    assert pin(root / name) == expected, name
prepared = {'commands': planned, 'source_inputs': inputs, 'source_proof': source,
    'authority': authority, 'backend_archive': str(archive_path.relative_to(root)),
    'object': str(object_path.relative_to(root)), 'observer': str(observer_path.relative_to(root)),
    'executable': str(host_path.relative_to(root)), 'frame_executable': str(frame_path.relative_to(root)),
    'header_dependency_inputs': header_names}
(work / 'prepared-build.json').write_text(json.dumps(prepared, indent=2) + '\n')
action = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'kind': 'diagnostic', 'expected_exit': 0,
    'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'argv': ['/usr/bin/python3', '-B', str(work / 'build_host.py')],
    'frozen_inputs': sorted(set(inputs) | {str((work / 'prepared-build.json').relative_to(root))}),
    'software_provider': False,
    'scope': 'Compile exactly1private retained_compute.cpp using the saved actualmatched7f2g++14flags; restore only1authenticated observer blob, replace only1archive member, relink original audited host and F frame with all other original inputs. M provider unchanged. No shader/ISPC/provider/observer/otherCPP builds, resources copied/restored or Sirius workloads. Source/build only; original runtime/scientific/release gates remain outstanding.'}
(work / 'build/actions.json').write_text(json.dumps({'host-build': action}, indent=2) + '\n')
print(json.dumps({'new_CPP_compilations': 1, 'new_links': 2, 'restored_observer_members': 1,
                  'frozen_inputs': len(inputs), 'compiled_or_executed': False}))
