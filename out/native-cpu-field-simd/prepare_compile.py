"""Seal six real compiler invocations; reuse the D/S compiler recipes exactly."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-field-simd'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


assert not os.environ.get('ISPC_ARGS')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (W / 'prepared-compile.json').exists()
assert not (W / 'compile/actions.json').exists()
assert not (B / 'objects').exists()
subprocess.run(['/usr/bin/python3', '-B', str(W / 'derive_executor.py')], check=True)
derived = doc(W / 'executor-derivation.json')
cpp = doc(R / 'out/native-cpu-static-program/batch-compile-result.json')
ispc = doc(R / 'out/native-cpu-zero-term/ispc-zero-compile-result.json')
assert ispc['all_objects_compiled'] and derived['field_width'] == 8
inputs = {**cpp['source_inputs'], **derived['source_inputs']}
objects = []
for slug, entry in derived['stages'].items():
    for kind in ('cpp_header', 'cpp_source', 'ispc_source'):
        item = entry[kind]
        inputs[item['path']] = {key: item[key] for key in ('bytes', 'sha256')}
    old = next(item for item in cpp['objects'] if Path(item['argv'][-3]).stem == slug)
    assert old['exit'] == 0 and old['argv'][-4] == '-c' and old['argv'][-2] == '-o'
    argv = list(old['argv'])
    include = argv.index(str(R / 'bin/linux-gcc/native-cpu-static-program/typed-batch'))
    assert argv[include-1] == '-isystem'
    argv[include] = str(B / 'cpp')
    argv[-3] = str(R / entry['cpp_source']['path'])
    output = 'bin/linux-gcc/native-cpu-field-simd/objects/' + slug + '-control.o'
    argv[-1] = str(R / output)
    inverse = list(argv)
    inverse[include], inverse[-3], inverse[-1] = old['argv'][include], old['argv'][-3], old['argv'][-1]
    assert inverse == old['argv']
    objects.append({'stage': slug, 'kind': 'cpp', 'output': output, 'argv': argv,
                    'original_argv': old['argv'], 'argv_inverse_exact': True})
    old = next(item for item in ispc['objects'] if item['stage'] == slug)
    assert old['exit'] == 0 and old['argv_executed']
    argv = list(old['argv'])
    argv[1] = str(R / entry['ispc_source']['path'])
    output = 'bin/linux-gcc/native-cpu-field-simd/objects/' + slug + '-layer.o'
    argv[-1] = str(R / output)
    inverse = list(argv); inverse[1], inverse[-1] = old['argv'][1], old['argv'][-1]
    assert inverse == old['argv']
    objects.append({'stage': slug, 'kind': 'ispc', 'output': output, 'argv': argv,
                    'original_argv': old['argv'], 'argv_inverse_exact': True})
tools = doc(R / 'out/native-cpu-host-representation/probe/preparation.json')['tool_dependencies']
for name, expected in tools.items():
    assert pin(Path(name)) == expected, name
    inputs[name] = expected
compiler = Path(next(item for item in objects if item['kind'] == 'ispc')['argv'][0])
assert pin(compiler)['sha256'] == '1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
inputs[str(compiler.relative_to(R))] = pin(compiler)
resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(compiler)], text=True, stderr=subprocess.STDOUT)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
for name in re.findall(r'(/[^\s()]+)', resolution):
    path = Path(name)
    if path.is_file():
        path = path.resolve(strict=True)
        inputs[str(path)] = pin(path)
for name in ('derive_executor.py', 'executor-derivation.json', 'prepare_compile.py', 'build.py',
             'compile/run_owned.py', 'compile/owned_guard.py'):
    inputs[str((W / name).relative_to(R))] = pin(W / name)
for name in ('out/native-cpu-static-program/batch-compile-result.json',
             'out/native-cpu-zero-term/ispc-zero-compile-result.json',
             'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             'out/native-cpu-row-ownership/cold-frame/owned_guard.py'):
    inputs[name] = pin(R / name)
for name, expected in inputs.items():
    assert pin(R / name) == expected, name
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
scope = 'Actual new all3 scalar-control bridges and all3 field-SIMD layer executors; original D/S compiler flags unchanged, all9380instructions and scalar row control retained. Source/build evidence only; no runtime, speed or scientific acceptance.'
prepared = {'revision': revision, 'inputs': inputs, 'objects': objects, 'scope': scope,
            'compiler_runtime_resolution': resolution,
            'transitive_header_limit': 'Original unchanged compiler recipes and direct headers sealed; no claim that every system transitive include is prefrozen.'}
(W / 'prepared-compile.json').write_text(json.dumps(prepared, indent=2) + '\n')
base = doc(R / 'out/native-cpu-zero-term/compile/actions.json')['ispc-zero']
actions = {}
for kind, seconds in (('cpp', 60), ('ispc', 300)):
    actions[kind] = {**base, 'revision': revision,
        'argv': ['/usr/bin/python3', '-B', str(W / 'build.py'), kind],
        'limits': {'seconds': seconds, 'rss_kib': 4194304, 'reserve_kib': 2097152},
        'frozen_inputs': sorted(set(inputs) | {str((W / 'prepared-compile.json').relative_to(R))}),
        'scope': scope}
(W / 'compile/actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'new_CPP_compilations_planned': 3, 'new_ISPC_compilations_planned': 3, 'executed': False}))
