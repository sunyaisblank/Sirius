"""Prepare one fixture compile/link, then one run; never execute either."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
P = R / 'out/native-cpu-upload-aware-staging'
B = R / 'bin/linux-gcc/native-cpu-upload-aware-boundary'
BACKEND = R / 'bin/linux-gcc/native-cpu-upload-aware-staging/matched/libsirius_backend.a'
LIMITS = {'seconds': 60, 'rss_kib': 4194304, 'reserve_kib': 2097152}


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2); stream.write('\n')


def owner_ok(path):
    d = json.loads(path.read_text())
    assert d['status'] == 'terminal' and d['passed'] and d['source_unchanged'] and d['whole_inputs_unchanged']
    assert not d['cleanup_errors'] and not d['remaining_owned_processes']
    for row in d['observed_births']:
        try:
            current = Path(f'/proc/{row["pid"]}/stat').read_text().rsplit(')', 1)[1].split()[19]
        except FileNotFoundError:
            current = None
        assert row['absent'] and current != row['saved_start_ticks']


mode, = sys.argv[1:]
assert mode in ('compile', 'run')
assert not any(os.environ.get(k) for k in ('CPATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH', 'GCC_EXEC_PREFIX', 'COMPILER_PATH'))
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
owner_ok(P / 'build/host-build/owner.json')
host = json.loads((P / 'host-build-result.json').read_text())
assert host['passed'] and host['source_inputs_unchanged']
assert R / host['backend_archive'] == BACKEND and pin(BACKEND) == host['backend_archive_identity']
inputs = dict(host['source_inputs'])
inputs[str(BACKEND.relative_to(R))] = pin(BACKEND)
for path in [P / 'host-build-result.json', P / 'build/host-build/owner.json', *W.glob('*.py'), W / 'boundary.cpp']:
    inputs[str(path.relative_to(R))] = pin(path)
if mode == 'compile':
    assert not (W / 'prepared.json').exists() and not (W / 'actions.json').exists()
    assert not (B / 'boundary.o').exists() and not (B / 'boundary').exists()
    source = next(c for c in host['commands'] if '-c' in c['argv'])
    link = next(c for c in host['commands'] if '-o' in c['argv'] and
                c['argv'][c['argv'].index('-o') + 1].endswith('/audited_backend_tests'))
    compile_argv = list(source['argv'])
    assert all(flag in compile_argv for flag in ('-O3', '-DNDEBUG', '-std=c++2c'))
    assert Path(compile_argv[0]).name.endswith('g++-14')
    compile_argv[compile_argv.index('-c') + 1] = str(W / 'boundary.cpp')
    compile_argv[compile_argv.index('-o') + 1] = str(B / 'boundary.o')
    link_argv, removed = [], []
    for word in link['argv']:
        if word.endswith('.o'):
            if not removed:
                link_argv.append(str(B / 'boundary.o'))
            removed.append(word)
        else:
            link_argv.append(word)
    assert removed and len([s for s in link_argv if s.endswith('.o')]) == 1
    assert link_argv.count(str(BACKEND)) == link['argv'].count(str(BACKEND)) > 0
    link_argv[link_argv.index('-o') + 1] = str(B / 'boundary')
    for command in (source, link):
        compiler = Path(command['argv'][0]).resolve(strict=True)
        inputs[str(compiler)] = pin(compiler)
    for word in link_argv:
        if word.endswith(('.a', '.so')) and not word.startswith('-'):
            path = Path(word)
            path = (path if path.is_absolute() else Path(link['cwd']) / path).resolve(strict=True)
            inputs[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    B.mkdir(parents=True, exist_ok=True)
    prepared = {'P_host_build_result': pin(P / 'host-build-result.json'), 'backend': str(BACKEND.relative_to(R)),
                'backend_identity': pin(BACKEND), 'source_inputs': inputs,
                'commands': [{'argv': compile_argv, 'cwd': source['cwd']}, {'argv': link_argv, 'cwd': link['cwd']}],
                'original_test_object_arguments_removed': removed, 'original_host_compile': source,
                'original_host_link': link, 'object': str((B / 'boundary.o').relative_to(R)),
                'executable': str((B / 'boundary').relative_to(R)),
                'scope': 'One extra diagnostic C++ fixture compile/link and fake wire-contract run. P actual host1CPP/two links and original benchmark observer remain separate and unchanged.'}
    write(W / 'prepared.json', prepared)
    actions = {}
else:
    owner_ok(W / 'compile/owner.json')
    built = json.loads((W / 'compile-result.json').read_text())
    assert built['passed'] and built['extra_fixture_CPP_compilations'] == built['extra_fixture_links'] == 1
    assert built['prepared'] == pin(W / 'prepared.json')
    inputs.update(json.loads((W / 'prepared.json').read_text())['source_inputs'])
    for key in ('object', 'executable'):
        inputs[built[key]['path']] = {k: built[key][k] for k in ('bytes', 'sha256')}
    inputs[str((W / 'compile-result.json').relative_to(R))] = pin(W / 'compile-result.json')
    inputs[str((W / 'compile/owner.json').relative_to(R))] = pin(W / 'compile/owner.json')
    actions = json.loads((W / 'actions.json').read_text())
    assert 'run' not in actions and not (W / 'run').exists()
inputs[str((W / 'prepared.json').relative_to(R))] = pin(W / 'prepared.json')
actions[mode] = {'revision': revision, 'kind': 'diagnostic', 'expected_exit': 0, 'software_provider': False,
                 'limits': LIMITS, 'frozen_inputs': sorted(inputs),
                 'argv': ['/usr/bin/python3', '-B', str(W / 'task.py'), mode, '{folder}'],
                 'scope': 'One owned upload-aware boundary ' + mode + '; actual P backend, fake returned-error/wire contract only; no math oracle, shader/device execution or adoption gate change.'}
if mode == 'compile':
    write(W / 'actions.json', actions)
else:
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared': mode, 'executed': False, 'limits': LIMITS}))
