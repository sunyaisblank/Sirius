"""Seal one dense C++ arithmetic build and its actual dependency closure."""
from pathlib import Path
import hashlib
import json
import re
import shlex
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-masked-arithmetic/probe'
assert sys.argv[1:] in ([], ['run'])


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
failure = R / 'out/native-cpu-dense-values-cohort/cold-frame'
owner = doc(failure / 'frame/owner.json')
assert owner['status'] == 'terminal' and not owner['passed'] and owner['stop_reason'] == 'diagnostic_timeout'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert not doc(failure / 'result.json')['accepted']
manifest = R / 'attestations/qualification' / owner['source']['revision'][:8] / 'native-cpu-dense-values-frame-timeout/MANIFEST.json'
assert not doc(manifest)['accepted']
assert subprocess.check_output(['git', 'show', '@{upstream}:' + str(manifest.relative_to(R))]) == manifest.read_bytes()
fp = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/drivers/llvmpipe/lp_state_cs.c'
assert pin(fp) == {'bytes': 109640, 'sha256': '34b1f093c09543c7f1c48873081d014849088b3839316dc6b3671d2b6add9cd1'}
prepared_path = W / 'prepared-cpp.json'
if not sys.argv[1:]:
    assert not prepared_path.exists() and not (W / 'cpp/actions.json').exists() and not B.exists()
    baseline = R / 'out/native-cpu-static-program/batch-compile-result.json'
    old = doc(baseline)['objects'][0]
    assert old['stage'] == 'Transport'
    flags = old['argv'][:old['argv'].index('-c')]
    assert flags.count('-march=x86-64') == 1 and '-std=c++2c' in flags
    assert '-frounding-math' in flags and '-ffp-contract=off' in flags
    entries = []
    inputs = {str(path.relative_to(R)): pin(path) for path in
              (baseline, manifest, fp, failure / 'frame/owner.json', failure / 'result.json',
               R / 'out/native-cpu-static-program/header-preparation.json', *W.glob('*.cpp'), *W.glob('*.py'))}
    for name in ('run_owned.py', 'owned_guard.py'):
        path = W / 'cpp' / name
        inputs[str(path.relative_to(R))] = pin(path)
    for slug in ('retained8',):
        argv = list(flags)
        if slug == 'retained8':
            argv[argv.index('-march=x86-64')] = '-march=x86-64-v3'
        source = W / (slug + '.cpp')
        output = B / (slug + '.o')
        dependency = subprocess.check_output([*argv, '-M', '-MT', 'dense-probe', str(source)], text=True)
        names = shlex.split(dependency.replace('\\\n', ' ').split(':', 1)[1])
        assert names
        for name in names:
            path = Path(name).resolve(strict=True)
            inputs[str(path)] = pin(path)
        entries.append({'name': slug, 'argv': [*argv, '-c', str(source), '-o', str(output)],
                        'object': str(output.relative_to(R)), 'dependency_words': len(names)})
    compiler = Path(flags[0]).resolve(strict=True)
    tools = {compiler, Path('/usr/bin/ld').resolve(), Path('/usr/bin/ldd').resolve()}
    tool_receipt = R / 'out/native-cpu-host-representation/probe/preparation.json'
    inputs[str(tool_receipt.relative_to(R))] = pin(tool_receipt)
    # Reuse discovery paths, while recording current bytes for this new build.
    # This does not repin or reinterpret the historical build receipt.
    for name in doc(tool_receipt)['tool_dependencies']:
        path = Path(name).resolve(strict=True)
        inputs[str(path)] = pin(path)
    for name in ('cc1plus', 'collect2'):
        path = Path(subprocess.check_output([str(compiler), '-print-prog-name=' + name], text=True).strip()).resolve(strict=True)
        assert path.is_file()
        tools.add(path)
    for path in list(tools):
        with path.open('rb') as stream:
            if stream.read(4) != b'\x7fELF':
                continue
        resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(path)], text=True, stderr=subprocess.STDOUT)
        assert 'undefined symbol' not in resolution and 'not found' not in resolution
        tools |= {Path(name).resolve() for name in re.findall(r'(/[^\s()]+)', resolution) if Path(name).is_file()}
    for path in tools:
        inputs[str(path)] = pin(path)
    driver_receipt = doc(R / 'out/native-cpu-dense-arithmetic/cpp-build-result.json')
    driver, = [row['object'] for row in driver_receipt['compilations'] if row['name'] == 'driver']
    assert pin(R / driver['path']) == {k:driver[k] for k in ('bytes','sha256')}
    inputs[driver['path']] = pin(R / driver['path'])
    inputs['out/native-cpu-dense-arithmetic/cpp-build-result.json'] = pin(R / 'out/native-cpu-dense-arithmetic/cpp-build-result.json')
    link = [str(compiler), str(R / driver['path']), str(B / 'retained8.o'), '-o', str(B / 'abi-check')]
    prepared = {'revision': revision, 'inputs': inputs, 'compilations': entries, 'link': link,
                'executable': str((B / 'abi-check').relative_to(R)),
                'scope': 'One new C++26 dense AVX2 helper compilation and link to the unchanged checked AC scalar-reference ABI driver object. Pure helpers are inlined within the two existing exports. Integer IEEE equality preserves NaN and signed-zero cases; operands in discarded arithmetic lanes are set to zero while selected ordered arithmetic remains unchanged. Existing exact D strict FP recipe retained; only helper uses guarded v3 to express AVX2 vectors. Driver stays baseline x86-64. Actual current compiler/runtime/header closure pinned; old D receipt supplies flags only, not current runtime identities. No ISPC/provider/host recompilation or workload in this build. Follow with one scoped mask/alias/refusal bitwise comparison, then exact source/code joins and fresh full38/C8/original frame before adoption. No independent math, renderer performance, science or release claim.'}
    prepared_path.write_text(json.dumps(prepared, indent=2) + '\n')
    inputs = set(inputs) | {str(prepared_path.relative_to(R))}
    kind, argv, label = 'build_only', ['/usr/bin/python3', '-B', str(W / 'build_cpp.py')], 'build'
else:
    prepared = doc(prepared_path)
    built = doc(W / 'cpp-build-result.json')
    old_owner = doc(W / 'cpp/build/owner.json')
    assert old_owner['status'] == 'terminal' and old_owner['passed']
    assert old_owner['source_unchanged'] and old_owner['whole_inputs_unchanged']
    assert not old_owner['cleanup_errors'] and not old_owner['remaining_owned_processes']
    assert built['passed'] and built['new_CPP_compilations'] == 1 and built['new_links'] == 1
    inputs = set(prepared['inputs']) | {str(prepared_path.relative_to(R))}
    inputs |= {str((W / name).relative_to(R)) for name in ('cpp-build-result.json', 'cpp/build/owner.json')}
    inputs |= set(built['runtime_libraries']) | {prepared['executable']}
    kind, argv, label = 'diagnostic', [str(R / prepared['executable'])], 'run'
action = {'revision': revision, 'kind': kind, 'expected_exit': 0,
          'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
          'argv': argv, 'frozen_inputs': sorted(inputs), 'software_provider': False,
          'scope': prepared['scope'] if label == 'build' else 'One new focused dense C++ arithmetic ABI comparison against original generated scalar helpers: whole scratch-bank words, inactive lanes, mixed masks, a/b/d aliasing, exact identities/signed zeros/subnormals/nonfinite/refusal boundaries. Scalar agreement is supporting finite compatibility, not independent math proof, full-program or performance acceptance.'}
actions = doc(W / 'cpp/actions.json') if (W / 'cpp/actions.json').exists() else {}
assert label not in actions
actions[label] = action
(W / 'cpp/actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared': label, 'executed': False, 'frozen_inputs': len(inputs)}))
