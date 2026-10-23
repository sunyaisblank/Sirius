"""Compile one selector; reuse every original observer and other GCC object."""
from pathlib import Path
import hashlib, json, shlex, shutil, subprocess
R = Path.cwd(); W = Path(__file__).resolve().parent
C = R / 'bin/linux-gcc'; B = C / 'native-cooperative-transport-intergroup/matched'
def sha(p): return hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()
modules = json.loads((W / 'compiler-result.json').read_text())['modules']
for record in modules.values():
    assert sha(R / record['module']) == record['sha256']
source = (W / 'retained_compute.cpp').read_text()
assert sha(W / 'retained_compute.cpp') == json.loads((W / 'source-scope.json').read_text())['host_selector']['sha256']
rows = json.loads((C / 'compile_commands.json').read_text())
row, = [x for x in rows if x['file'].endswith('/backend/retained_compute.cpp')]
args = shlex.split(row['command'])
args[args.index('-o') + 1] = str(B / 'retained_compute.cpp.o')
args[args.index('-c') + 1] = str(W / 'retained_compute.cpp')
args[1:1] = ['-I' + str(W / 'headers')]
commands = [{'argv': args, 'cwd': row['directory']}]
subprocess.run(args, cwd=row['directory'], check=True)
original = C / 'src/sirius/backend/libsirius_backend.a'
candidate = B / 'libsirius_backend.a'; before = sha(original)
shutil.copyfile(original, candidate)
subprocess.run(['/usr/bin/ar', 'rDs', str(candidate), str(B / 'retained_compute.cpp.o')], check=True)
members = subprocess.check_output(['/usr/bin/ar', 't', str(original)], text=True).splitlines()
assert members == subprocess.check_output(['/usr/bin/ar', 't', str(candidate)], text=True).splitlines()
different = []
for member in members:
    old = subprocess.check_output(['/usr/bin/ar', 'p', str(original), member])
    new = subprocess.check_output(['/usr/bin/ar', 'p', str(candidate), member])
    if old != new: different.append(member)
assert different == ['retained_compute.cpp.o']
lines = subprocess.check_output(['ninja', '-C', str(C), '-t', 'commands',
                                 'tests/backend/sirius_backend_tests'], text=True).splitlines()
line, = [x for x in lines if ' -o tests/backend/sirius_backend_tests ' in x]
args = shlex.split(line); start = args.index('/usr/bin/g++-14'); end = args.index('&&', start)
args = args[start:end]
old = 'tests/backend/CMakeFiles/sirius_backend_tests.dir/retained_compute_test.cpp.o'
assert args.count(old) == 1
args[args.index(old)] = str(R / 'bin/linux-gcc/native-original-input-ownership/matched/audited-observer.o')
args[args.index('-o') + 1] = str(B / 'audited_backend_tests')
assert args.count('src/sirius/backend/libsirius_backend.a') == 2
args = [str(candidate) if x == 'src/sirius/backend/libsirius_backend.a' else x for x in args]
commands.append({'argv': args, 'cwd': str(C)})
subprocess.run(args, cwd=C, check=True)
assert sha(original) == before
for record in modules.values():
    data = (R / record['module']).read_bytes()
    assert (B / 'retained_compute.cpp.o').read_bytes().count(data) == 1
    assert (B / 'audited_backend_tests').read_bytes().count(data) == 1
(W / 'matched-result.json').write_text(json.dumps({
    'passed': True, 'commands': commands, 'new_cpp_compilations': 1,
    'observer_compilations': 0, 'changed_archive_members': different,
    'original_backend_archive_unchanged': True,
    'original_public_source': '0109ab462edade524667d20303507aa09adcbaf6',
    'same_actual_executable_for_both_provider_arms': True,
    'executable': str((B / 'audited_backend_tests').relative_to(R)),
    'executable_sha256': sha(B / 'audited_backend_tests'),
    'actual_three_module_spans_present_once_in_selector_object_and_final_ELF': True,
    'scope': 'Three software-CPU-only native32 selectors only; all35 original arrays and every wide/FMA/portable/owner route preserved. No numerical, integrated, benefit, frame or release acceptance.'}, indent=2) + '\n')
print(json.dumps({'passed': True, 'new_cpp_units': 1, 'observer_compilations': 0,
                  'changed_archive_members': different}))
