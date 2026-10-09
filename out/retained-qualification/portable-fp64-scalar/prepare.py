#!/usr/bin/env python3
"""Offline preparation; never enumerates or dispatches a compute device."""
import collections
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(args):
    completed = subprocess.run(list(map(str, args)), cwd=ROOT, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(f"command failed {args}:\n{completed.stdout}\n{completed.stderr}")
    return completed.stdout

fixture = ROOT / 'bin/linux-gcc/tests/backend/portable_binary32_reference.bin'
reference = ROOT / 'tests/support/portable_binary32/reference.py'
canonical = ROOT / 'src/sirius/kernels/portable_binary32.h'
assert digest(canonical) == '2f835dc2cbda7a62425977b31dc89155bc0393284d5ce3ace90f42c71bd076f2'
assert digest(fixture) == '3fb1433f9cd14fbc9fd2b99a38a4893a4a571150df10f2fbe58172318fc7691e'
blob = fixture.read_bytes()
assert struct.unpack_from('<III', blob) == (0x32334250, 1, 174268)
assert len(blob) == 12 + 174268 * 24
inputs, expected, counts = [174268], [], collections.Counter()
for op, a, b, high, low, valid in struct.iter_unpack('<6I', blob[12:]):
    assert 0 <= op <= 13 and valid <= 1
    assert op == 13 or (low == 0 and valid == 1)
    inputs.extend((op, a, b))
    expected.extend((high, low, valid))
    counts[op] += 1
assert counts[13] == 53524 and sum(counts.values()) == 174268
assert sum(n for op, n in counts.items() if op != 13) == 120744
for name, words in [('shader-input.bin', inputs), ('shader-expected.bin', expected)]:
    (HERE / name).write_bytes(struct.pack('<' + str(len(words)) + 'I', *words))
assert digest(HERE/'shader-input.bin') == '1687719c972f8eac46b22b7f94513893413561c0dd1c612ff62f8297bef6ceac'
assert digest(HERE/'shader-expected.bin') == '293b4c795a20cd75556697564af0e653788c734f032326f812cc902f85ff25e6'

compile_shader = ['/opt/slang/bin/slangc', HERE/'raw_word_probe.slang', '-I', ROOT/'src',
                  '-I', HERE, '-O3', '-target', 'spirv', '-profile', 'spirv_1_5',
                  '-entry', 'ComputeMain', '-stage', 'compute', '-o', HERE/'probe.raw.spv']
run(compile_shader)
run(['spirv-dis', HERE/'probe.raw.spv', '-o', HERE/'probe.raw.spvasm'])
assembly = (HERE/'probe.raw.spvasm').read_text()
entry = re.search(r'OpEntryPoint GLCompute (%\S+)', assembly)[1]
assert 'OpCapability RoundingModeRTE' not in assembly
assembly = assembly.replace('OpCapability Shader', 'OpCapability Shader\n               OpCapability RoundingModeRTE', 1)
local_size = re.search(r'^.*OpExecutionMode ' + re.escape(entry) + r' LocalSize.*$', assembly, re.M)[0]
assert local_size.strip() == f'OpExecutionMode {entry} LocalSize 64 1 1'
assembly = assembly.replace(local_size, local_size + '\n               OpExecutionMode ' + entry + ' RoundingModeRTE 64', 1)
(HERE/'probe.patched.spvasm').write_text(assembly)
run(['spirv-as', '--target-env', 'spv1.5', HERE/'probe.patched.spvasm', '-o', HERE/'probe.spv'])
run(['spirv-val', '--target-env', 'vulkan1.2', HERE/'probe.spv'])
# Fresh disassembly of the actual final module; do not inspect a stale neighbor.
run(['spirv-dis', HERE/'probe.spv', '-o', HERE/'probe.spvasm'])
assembly = (HERE/'probe.spvasm').read_text()
capabilities = re.findall(r'OpCapability (\S+)', assembly)
assert set(capabilities) == {'Shader', 'Float64', 'RoundingModeRTE'} and len(capabilities) == 3
assert set(re.findall(r'OpTypeInt (\d+) [01]', assembly)) == {'32'}
assert re.findall(r'OpTypeFloat (\d+)', assembly) == ['64']
assert ' Fma ' not in assembly and 'OpFConvert' not in assembly
float_operations = re.findall(r'(%\S+) = (OpF(?!unction)\w+) (\S+)', assembly)
double_type = re.search(r'(%\S+) = OpTypeFloat 64', assembly)[1]
assert float_operations and {op for _, op, _ in float_operations} == {'OpFAdd', 'OpFMul'}
assert all(t == double_type for _, _, t in float_operations)
decorated = set(re.findall(r'OpDecorate (%\S+) NoContraction', assembly))
assert all(result in decorated for result, _, _ in float_operations)
entry = re.search(r'OpEntryPoint GLCompute (%\S+)', assembly)[1]
modes = re.findall(r'OpExecutionMode\S* (%\S+) (.*)', assembly)
assert set(modes) == {(entry, 'LocalSize 64 1 1'), (entry, 'RoundingModeRTE 64')}
assert set(re.findall(r'OpDecorate %\S+ DescriptorSet (\d+)', assembly)) == {'0'}
assert sorted(re.findall(r'OpDecorate %\S+ Binding (\d+)', assembly)) == ['0', '1']
assert set(re.findall(r'OpDecorate %\S+ ArrayStride (\d+)', assembly)) == {'4'}

common = ['/usr/bin/g++-14', '-std=c++2c', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-Werror',
          '-fno-fast-math', '-ffp-contract=off', '-I', ROOT/'src', '-I', HERE]
compile_host = common + ['-frounding-math', HERE/'host_probe.cpp', '-o', HERE/'host_probe']
run(compile_host)
libs = [ROOT/'bin/linux-gcc/src/sirius/backend/libsirius_backend.a',
        ROOT/'bin/linux-gcc/src/sirius/core/libsirius_core.a',
        ROOT/'bin/linux-gcc/src/sirius/libsirius_base.a']
compile_gpu = common + ['-DSIRIUS_CONTRACT_MODE=2', '-DSIRIUS_HAS_RETAINED_COMPUTE=1',
                       HERE/'gpu_probe.cpp', '-o', HERE/'gpu_probe', *libs, '-lvulkan', '-ldl', '-pthread']
link_library_hashes = {str(p.relative_to(ROOT)): digest(p) for p in libs}
run(compile_gpu)
assert all(digest(ROOT/name) == value for name, value in link_library_hashes.items()), 'link library changed during diagnostic link'
host_result = json.loads(run([HERE/'host_probe', HERE/'shader-input.bin', HERE/'shader-expected.bin', HERE/'host-actual.bin']))
assert host_result == {'scope': 'isolated host raw-word prototype', 'cases': 174268,
                       'observed_words': 522804, 'mismatches': 0}
assert digest(HERE/'host-actual.bin') == digest(HERE/'shader-expected.bin')
preflight = json.loads(run([HERE/'gpu_probe', HERE/'probe.spv', HERE/'shader-input.bin', HERE/'shader-expected.bin', HERE/'gpu-actual.bin']))
assert preflight['device_calls'] == preflight['dispatches'] == 0
report = {
    'scope': 'isolated finite raw-word prototype; no production integration or retained-stage qualification',
    'source_head': run(['git', 'rev-parse', 'HEAD']).strip(),
    'branch': run(['git', 'branch', '--show-current']).strip(),
    'tracked_working_tree': run(['git', 'status', '--short', '--untracked-files=no']),
    'corpus': {'path': str(fixture.relative_to(ROOT)), 'sha256': digest(fixture),
               'reference_source': str(reference.relative_to(ROOT)), 'reference_source_sha256': digest(reference),
               'extraction': 'strict parse and exact word split; no expected arithmetic recomputed',
               'records': 174268, 'primitive_records': 120744, 'residual_records': 53524,
               'operation_counts': dict(sorted(counts.items()))},
    'artifact_guards': {'validator': 'spirv-val --target-env vulkan1.2', 'capabilities': capabilities,
                       'integer_widths': [32], 'float_widths': [64], 'rounding': 'RoundingModeRTE64 only',
                       'native_double_operations': float_operations, 'all_NoContraction': True,
                       'no_native_float32': True, 'no_Int64': True, 'array_stride': 4,
                       'descriptor_bindings': [0, 1], 'workgroup': [64, 1, 1]},
    'commands': {'shader': list(map(str, compile_shader)), 'host': list(map(str, compile_host)),
                 'gpu_harness': list(map(str, compile_gpu))},
    'tool_versions': {'slang': run(['/opt/slang/bin/slangc', '-version']).strip(),
                      'g++': run(['/usr/bin/g++-14', '--version']).splitlines()[0],
                      'spirv-val': run(['spirv-val', '--version']).strip()},
    'host_result': host_result, 'gpu_preflight': preflight,
    'linked_libraries_sha256': link_library_hashes,
    'bound_sources': {str(p.relative_to(ROOT)): digest(p) for p in
                     [canonical, reference, fixture, *libs, *[HERE/x for x in
                       ['fp64_scalar.h', 'raw_word_probe.slang', 'word_io.h', 'host_probe.cpp', 'gpu_probe.cpp',
                        'prepare.py', 'probe.spv', 'probe.spvasm', 'shader-input.bin', 'shader-expected.bin',
                        'host-actual.bin', 'host_probe', 'gpu_probe']]]},
}
(HERE/'preparation-report.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'prepared': True, 'host': host_result, 'module_sha256': digest(HERE/'probe.spv'),
                  'assisted_operation_records': sum(counts[x] for x in [0, 1, 2]), 'device_calls': 0}))
