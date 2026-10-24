"""Prepare compile, then run only after terminal owned compile; execute neither."""
from pathlib import Path
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
EXE = R / 'bin/linux-gcc/native-cpu-host-representation/spacing-probe'
LIMITS = {'seconds': 60, 'rss_kib': 4194304, 'reserve_kib': 2097152}


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2); stream.write('\n')


def check_owner(path):
    owner = json.loads(path.read_text())
    assert owner['status'] == 'terminal' and owner['passed']
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not any(os.environ.get(key) for key in
               ('CPATH', 'CPLUS_INCLUDE_PATH', 'C_INCLUDE_PATH', 'LIBRARY_PATH', 'GCC_EXEC_PREFIX',
                'COMPILER_PATH', 'GCC_SPECS', 'LD_PRELOAD', 'LD_AUDIT', 'LD_LIBRARY_PATH'))
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
mode, = sys.argv[1:]
assert mode in ('compile', 'run')
generation = json.loads((W / 'generation.json').read_text())
assert pin(W / 'probe.cpp') == generation['generated_source']
for key in ('source', 'patch'):
    assert pin(R / generation[key]) == generation[key + '_identity']
inputs = {str(path.relative_to(R)) for path in W.glob('*.py')}
inputs |= {'src/sirius/backend/retained_compute.cpp', str((W.parent / 'host.patch').relative_to(R)),
           str((W.parent / 'proof.json').relative_to(R)), str((W / 'probe.cpp').relative_to(R)),
           str((W / 'generation.json').relative_to(R)), 'bin/linux-gcc/compile_commands.json',
           'out/native-cpu-row-ownership/cold-frame/run_owned.py',
           'out/native-cpu-row-ownership/cold-frame/owned_guard.py'}
if mode == 'compile':
    assert not (W / 'actions.json').exists() and not EXE.exists()
    command = next(row['command'] for row in json.loads((R / 'bin/linux-gcc/compile_commands.json').read_text())
                   if row['file'] == str(R / 'src/sirius/backend/retained_compute.cpp'))
    tokens = shlex.split(command); compiler = Path(tokens[0]).resolve(strict=True)
    assert Path(tokens[0]).name == 'g++-14'
    assert all(flag in tokens for flag in ('-O3', '-DNDEBUG', '-std=c++2c'))
    assert not any(flag.startswith(('-ffast-math', '-Ofast', '-funsafe-math', '-ffinite-math', '-fno-builtin'))
                   for flag in tokens)
    def query(flag):
        return subprocess.check_output([str(compiler), flag], text=True, timeout=10).strip()
    def resolve_tool(value):
        path = Path(value)
        return (path if path.is_absolute() else Path(shutil.which(value))).resolve(strict=True)
    tools = {str(compiler): pin(compiler), str(Path(sys.executable).resolve()): pin(Path(sys.executable).resolve())}
    for name in ('cc1plus', 'collect2', 'as', 'ld'):
        path = resolve_tool(query('-print-prog-name=' + name)); tools[str(path)] = pin(path)
    for name in ('libstdc++.so', 'libgcc_s.so.1', 'libm.so.6', 'libc.so.6', 'ld-linux-x86-64.so.2'):
        value = query('-print-file-name=' + name)
        path = Path(value) if value != name else Path('/lib/x86_64-linux-gnu') / name
        path = path.resolve(strict=True); tools[str(path)] = pin(path)
    ldd = Path(shutil.which('ldd')).resolve(strict=True); tools[str(ldd)] = pin(ldd)
    direct_headers = {}
    for name in ('algorithm', 'array', 'bit', 'cerrno', 'cfenv', 'cmath', 'cstdint', 'cstdio',
                 'cstring', 'limits', 'stdexcept', 'vector'):
        path = Path('/usr/include/c++/14') / name; direct_headers[str(path)] = pin(path)
    include = Path(query('-print-file-name=include'))
    for path in (include / 'immintrin.h', include / 'xmmintrin.h', include / 'fxsrintrin.h',
                 Path('/usr/include/fenv.h'), Path('/usr/include/x86_64-linux-gnu/bits/fenv.h')):
        direct_headers[str(path)] = pin(path)
    argv, position = [str(compiler)], 1
    while position < len(tokens):
        word = tokens[position]; position += 1
        if word == '-o':
            position += 1
        elif word == '-c':
            pass
        elif word == str(R / 'src/sirius/backend/retained_compute.cpp'):
            argv.append(str(W / 'probe.cpp'))
        else:
            argv.append(word)
    argv += ['-o', str(EXE), '-MD', '-MF', str(W / 'compile/dependencies.d'), '-Wl,-t', '-lm']
    EXE.parent.mkdir(exist_ok=True)
    preparation = {'revision': revision, 'genuine_host_compile_command': command, 'argv': argv,
                   'compiler_version': query('-dumpfullversion'), 'tool_dependencies': tools,
                   'direct_headers': direct_headers, 'ldd': str(ldd), 'executable': str(EXE.relative_to(R)),
                   'limits': LIMITS, 'scope': 'One actual-lambda numeric/errno/fenv/MXCSR boundary probe; no performance, physics or qualification claim.',
                   'expected_unique_inputs': 2572, 'expected_environments': 32, 'expected_comparisons': 82304,
                   'transitive_header_limit': 'Single compile emits -MD; transitive header identities are recorded afterward, not falsely claimed pre-frozen.'}
    write(W / 'preparation.json', preparation)
    inputs |= set(tools) | set(direct_headers) | {str((W / 'preparation.json').relative_to(R))}
    actions = {}
else:
    check_owner(W / 'compile/owner.json')
    receipt = json.loads((W / 'compile-result.json').read_text())
    assert receipt['passed'] and pin(EXE) == receipt['executable']
    inputs |= set(receipt['headers']) | set(receipt['runtime_libraries']) | set(receipt['link_inputs'])
    inputs |= set(json.loads((W / 'preparation.json').read_text())['tool_dependencies'])
    inputs |= {str(EXE.relative_to(R)), str((W / 'preparation.json').relative_to(R)),
               str((W / 'compile-result.json').relative_to(R))}
    actions = json.loads((W / 'actions.json').read_text()); assert 'run' not in actions
    assert not (W / 'run').exists()
actions[mode] = {'revision': revision, 'kind': 'build', 'expected_exit': 0, 'software_provider': False,
                 'limits': LIMITS, 'scope': 'One owned spacing ' + mode + '; fixed boundary/side-effect matrix, not a performance or physics study.',
                 'argv': [str(Path(sys.executable).resolve()), '-B', str(W / 'task.py'), mode, '{folder}'],
                 'frozen_inputs': sorted(inputs)}
if mode == 'compile':
    write(W / 'actions.json', actions)
else:
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared': mode, 'executed': False, 'limits': LIMITS}))
