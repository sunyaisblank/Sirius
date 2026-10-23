"""Seal one selector compilation/link with all original observer and library inputs."""
from pathlib import Path
import hashlib, json, re, shlex, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent; C = R / 'bin/linux-gcc'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
scope = doc(W / 'source-scope.json'); assert scope['source_review_passed']
assert pin(W / 'retained_compute.cpp') == scope['host_selector']
assert pin(W / 'headers/retained_kernels.h') == scope['header']
inputs = {str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in (C / 'src/sirius/backend/retained').iterdir() if p.is_file()}
for record in doc(W / 'compiler-result.json')['modules'].values():
    assert pin(R / record['module']) == {k: record[k] for k in ['bytes', 'sha256']}
    inputs.add(record['module'])
inputs |= {'bin/linux-gcc/CMakeCache.txt', 'bin/linux-gcc/compile_commands.json', 'bin/linux-gcc/build.ninja',
    'bin/linux-gcc/native-original-input-ownership/matched/audited-observer.o',
    'out/native-cpu-row-ownership/cold-frame/run_owned.py',
    'out/native-cpu-row-ownership/cold-frame/owned_guard.py', '/usr/bin/python3', '/etc/ld.so.cache'}
observer = R / 'bin/linux-gcc/native-original-input-ownership/matched/audited-observer.o'
assert pin(observer) == {'bytes': 9376200, 'sha256': '4d9f8dd3364c08a58fb8ed56540331fecd1d94a66f221024c6c91f602cd72355'}
lines = subprocess.check_output(['ninja', '-C', str(C), '-t', 'commands',
    'tests/backend/sirius_backend_tests'], text=True).splitlines()
line, = [x for x in lines if ' -o tests/backend/sirius_backend_tests ' in x]
for word in shlex.split(line):
    p = Path(word); p = p if p.is_absolute() else C / p
    if p.is_file(): inputs.add(str(p.relative_to(R)) if p.is_relative_to(R) else str(p.resolve()))
frontend = subprocess.check_output(['/usr/bin/g++-14', '-print-prog-name=cc1plus'], text=True).strip()
for tool in ['/usr/bin/g++-14', '/usr/bin/ar', '/usr/bin/ninja', '/usr/bin/as', '/usr/bin/ld', frontend]:
    p = Path(tool); assert p.is_file(); inputs.add(str(p.resolve()))
    result = subprocess.run(['ldd', str(p)], capture_output=True, text=True)
    inputs |= {str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)', result.stdout) if Path(x).is_file()}
for name in inputs: assert (R / name).is_file(), name
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
action = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'kind': 'diagnostic', 'expected_exit': 0,
    'limits': {'seconds': 180, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'argv': ['/usr/bin/python3', '-B', str(W / 'build_matched.py')],
    'frozen_inputs': sorted(inputs), 'software_provider': False,
    'scope': 'One private selector object compilation and original observer executable link, composing three immutable SPVs. Original observer and every other archive member/object retained; actual object/ELF module data joins required. No shader/provider compilation, runtime/numerical/benefit/frame/scientific/release acceptance.'}
target = W / 'actions.json'; assert not target.exists()
target.write_text(json.dumps({'matched': action}, indent=2) + '\n')
print(json.dumps({'frozen_inputs': len(inputs), 'new_CPP_units_planned': 1,
                  'new_shader_provider_observer_compilations': 0}))
