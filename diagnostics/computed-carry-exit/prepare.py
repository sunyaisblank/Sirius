"""Focused issue101 recipe: pin, compile four probes, or verify finished owners.

No mode invokes a provider. Root runs each explicit provider action separately.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
BUILD = ROOT / 'bin/linux-gcc/computed-carry-exit-prerequisites'
PAIRS = [(p, m) for p in ('stock', 'candidate') for m in ('integer', 'normal-transport')]


def record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}


def write(name, data):
    (WORK / name).write_text(json.dumps(data, indent=2) + '\n')


def frozen():
    for name, expected in json.loads((WORK / 'source-freeze.json').read_text()).items():
        assert record(name) == expected, name


def pin():
    frozen()
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'])
    assert not BUILD.exists() and not (WORK / 'operator-context.json').exists()
    actions = json.loads((WORK / 'actions.json').read_text())
    assert all(not (WORK / name).exists() for name in actions)
    assert all(a['revision'] == 'PENDING_OPERATOR_CLEAN_CONTEXT' for a in actions.values())
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    for action in actions.values():
        action['revision'] = revision
    write('actions.json', actions)
    write('operator-context.json', {'revision': revision, 'status': '',
        'actions': record(WORK / 'actions.json'), 'source_freeze': record(WORK / 'source-freeze.json')})
    print('Pinned clean operator context', revision, flush=True)


def compile_probes():
    frozen()
    context = json.loads((WORK / 'operator-context.json').read_text())
    assert context['actions'] == record(WORK / 'actions.json')
    assert context['source_freeze'] == record(WORK / 'source-freeze.json')
    assert not BUILD.exists()
    BUILD.mkdir()
    cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
    tool = lambda n: re.search(r'^' + n + r':FILEPATH=(.*)$', cache, re.M)[1]
    compiler, optimizer, disassembler, assembler, validator = [tool(n) for n in
        ('SIRIUS_SLANGC', 'SIRIUS_SPIRV_OPT', 'SIRIUS_SPIRV_DIS', 'SIRIUS_SPIRV_AS', 'SIRIUS_SPIRV_VAL')]
    spec = importlib.util.spec_from_file_location('issue101_original_policy', ROOT / 'scripts/build-retained-kernels.py')
    policy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(policy)
    modules = {}
    for producer, mode in PAIRS:
        stem = producer + '-' + mode
        raw, optimized = BUILD / (stem + '.raw.spv'), BUILD / (stem + '.optimized.spv')
        assembly, output = BUILD / (stem + '.spvasm'), BUILD / (stem + '.spv')
        definitions = ['-DSIRIUS_RETAINED_PORTABLE=1']
        if mode == 'normal-transport':
            definitions += ['-DSIRIUS_RETAINED_NORMAL_SUM32=1', '-DSIRIUS_RETAINED_PARALLEL_TRANSPORT=1']
        subprocess.run([compiler, str(WORK / 'probe.slang'), *definitions, '-I', str(WORK / producer),
            '-O0', '-target', 'spirv', '-profile', 'spirv_1_5', '-entry', 'ComputeMain',
            '-stage', 'compute', '-o', str(raw)], check=True)
        subprocess.run([optimizer, '--target-env=vulkan1.2', '--ssa-rewrite',
            '--eliminate-dead-code-aggressive', '--preserve-bindings', '--preserve-interface',
            str(raw), '-o', str(optimized)], check=True)
        subprocess.run([disassembler, str(optimized), '-o', str(assembly)], check=True)
        text = assembly.read_text()
        if mode == 'normal-transport':
            text = policy.portable_normal_sum_controls(text)
            assembly.write_text(text)
        else:
            assert re.findall(r'OpCapability (\S+)', text) == ['Shader']
            assert set(re.findall(r'OpTypeInt (\d+) [01]', text)) == {'32'}
            assert 'OpTypeFloat' not in text
        assert ' LocalSize 64 1 1' in text
        assert len(re.findall(r'OpDecorate %\S+ Binding [01]$', text, re.M)) == 2
        assert 'RTAccumulateTerm' in text
        subprocess.run([assembler, '--target-env', 'vulkan1.2', str(assembly), '-o', str(output)], check=True)
        subprocess.run([validator, '--target-env', 'vulkan1.2', str(output)], check=True)
        modules[stem] = record(output)
        print(stem, 'complete and validated', flush=True)
    command = ['/usr/bin/g++-14', '-DSIRIUS_CONTRACT_MODE=2', '-DSIRIUS_HAS_RETAINED_COMPUTE=1',
        '-DSIRIUS_RELEASE_RESOURCE_LOCKED=1', '-I' + str(ROOT / 'src'), '-O3', '-DNDEBUG', '-std=c++2c',
        '-Wall', '-Wextra', '-Wpedantic', '-Wno-missing-field-initializers', '-Werror',
        '-MMD', '-MF', str(BUILD / 'gpu.d'), '-c', str(WORK / 'gpu.cpp'), '-o', str(BUILD / 'gpu.o')]
    subprocess.run(command, check=True)
    subprocess.run(['/usr/bin/g++-14', str(BUILD / 'gpu.o'),
        str(ROOT / 'bin/linux-gcc/src/sirius/backend/CMakeFiles/sirius_backend.dir/vulkan/vulkan_device.cpp.o'),
        str(ROOT / 'bin/linux-gcc/src/sirius/libsirius_base.a'), '-lvulkan', '-ldl', '-pthread',
        '-o', str(BUILD / 'gpu')], check=True)
    frozen()
    write('compile-result.json', {'modules': modules, 'gpu': record(BUILD / 'gpu'),
        'scope': 'Four complete actual-function probes and one harness compiled/validated; no provider or arithmetic execution.'})


def verify_finished():
    frozen()
    context = json.loads((WORK / 'operator-context.json').read_text())
    expected = (WORK / 'expected.bin').read_bytes()
    reports = {}
    for name in ['compile-focused'] + ['provider-' + p + '-' + m for p, m in PAIRS]:
        folder = WORK / name
        owner = json.loads((folder / 'owner.json').read_text())
        assert owner['passed'] is True and owner['child_exit'] == 0 and owner['status'] == 'terminal'
        assert owner['source'] == {'revision': context['revision'], 'status': ''}
        assert owner['stop_reason'] is None and owner['cleanup_errors'] == []
        assert owner['remaining_owned_processes'] == []
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert all(b['absent'] for b in owner['observed_births'])
        assert json.loads((folder / 'inputs-before.json').read_text()) == json.loads((folder / 'inputs-after.json').read_text())
        if name.startswith('provider-'):
            assert owner['provider_mappings']
            actual = (folder / 'actual.bin').read_bytes()
            assert actual == expected
            summary = [json.loads(line) for line in (folder / 'stdout.log').read_text().splitlines()
                       if line.startswith('{') and 'completed_dispatches' in line]
            assert len(summary) == 1 and summary[0]['completed_dispatches'] == 1
            assert summary[0]['mismatches'] == summary[0]['untouched_complement_words'] == summary[0]['malformed_words'] == 0
            reports[name] = {'actual': record(folder / 'actual.bin'), 'summary': summary[0]}
    write('focused-result.json', {'passed': True, 'providers': reports,
        'fixture': record(WORK / 'fixture.json'), 'scope': 'Focused complete five-word stock/candidate comparison against independent Fraction expectations in two modes only. Original174268 corpus, eight integrated controls, all products/consumers and benefit gates remain deferred and required.'})


if __name__ == '__main__':
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024**2, 64 * 1024**2))
    {'pin-context': pin, 'compile': compile_probes, 'verify': verify_finished}[sys.argv[1]]()
