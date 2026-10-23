"""Prepare the three typed CPU adapters and complete-table admission; no compiler."""
from pathlib import Path
import hashlib
import json
import re
import shutil

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


batch = json.loads((W / 'batch-generation.json').read_text())
authority = json.loads((R / 'out/native-cpu-row-ownership/compiler-result.json').read_text())
table_path = R / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert pin(table_path) == authority['original_inputs'][str(table_path.relative_to(R))]
table_source = table_path.read_text()
policies = {
    'Transport': ('transport', 15967, 230, 4699, 9696,
                  {'retainedScratch_0': 8, 'retainedStatus_0': 9188,
                   'retainedCoefficients_0': 9192}),
    'Endpoint': ('endpoint', 21143, 225, 3829, 12280,
                 {'nativeProjection_0': 8, 'retainedStatus_0': 36,
                  'retainedScratch_0': 40}),
    'DopriPhase': ('dopri_phase', 10689, 1810, 4149, 12592,
                   {'retainedStatus_0': 8, 'retainedScratch_0': 12}),
}
compiler = str(Path(shutil.which('g++-14')).resolve(strict=True))
inputs = [table_path, W / 'batch-generation.json', W / 'prepare_batch.py',
          W / 'native_packet.hpp', W / 'native_bridge.h']
prepared = []
B.mkdir(parents=True, exist_ok=False)
for stage, (slug, words, input_row, output_row, context_size, offsets) in policies.items():
    item = batch['modules'][stage]
    original = R / item['source']
    assert pin(original) == {k: item[k] for k in ('bytes', 'sha256')}
    reflection = R / item['reflection']
    reflected = json.loads(reflection.read_text())
    assert reflected['entryPoints'][0]['name'] == 'ComputeMain'
    assert reflected['entryPoints'][0]['stage'] == 'compute'
    assert reflected['entryPoints'][0]['threadGroupSize'] == [64, 1, 1]
    assert len(reflected['parameters']) == 2
    for index, name in enumerate(('inputs', 'outputs')):
        parameter = reflected['parameters'][index]
        assert parameter['name'] == name
        assert parameter['binding'] == {
            'kind': 'uniform', 'offset': index * 16, 'size': 16, 'elementStride': 0}
        assert parameter['type']['baseShape'] == 'structuredBuffer'
        assert parameter['type']['resultType'] == {'kind': 'scalar', 'scalarType': 'uint32'}
    assert reflected['parameters'][1]['type']['access'] == 'readWrite'
    raw = original.read_text()
    prelude = '#include "/opt/slang/include/slang-cpp-prelude.h"\n'
    assert raw.startswith(prelude)
    assert raw.count('precise float ') == 4
    fixed = raw[len(prelude):].replace('precise float ', 'float ')
    fixed, count = re.subn(r'\bKernelContext_0 (kernelContext_\d+);',
                           r'KernelContext_0 \1{};', fixed)
    assert count == 1
    entry_names = ('_ComputeMain', 'ComputeMain_Thread', 'ComputeMain_Group', 'ComputeMain')
    for name in entry_names:
        assert re.search(r'\bvoid ' + name + r'\(', fixed), name
        fixed = re.sub(r'\b' + name + r'\b', 'sirius_' + slug + name, fixed)
    namespace = 'sirius_generated_' + slug
    header = B / (slug + '.native.hpp')
    header.write_text(prelude + '\nnamespace ' + namespace + ' {\n' + fixed + '\n}\n')
    match = re.search(r'inline constexpr std::array<std::uint32_t, ' + str(words) +
                      r'> k' + stage + r'Program\{\{(.*?)\}\};', table_source, re.S)
    assert match
    values = [int(value) for value in re.findall(r'(\d+)u', match[1])]
    assert len(values) == words and all(0 <= value < 2**32 for value in values)
    assert not re.sub(r'[\d\s,u]', '', match[1])
    canonical = B / (slug + '.program.hpp')
    canonical.write_text('#pragma once\n#include <cstdint>\nnamespace sirius_cpu_bridge {\n'
                         'inline constexpr std::uint32_t ' + slug + '_program[] = {\n' +
                         match[1] + '\n};\n}\n')
    assertions = '\n'.join('static_assert(offsetof(g::KernelContext_0, ' + field +
                            ') == ' + str(offset) + ');' for field, offset in offsets.items())
    unit = B / (slug + '.cpp')
    unit.write_text('''#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include "native_bridge.h"
#include "native_packet.hpp"
#include <SLUG.native.hpp>
#include <SLUG.program.hpp>
namespace g = NAMESPACE;
static_assert(sizeof(void*) == 8 && sizeof(float) == 4);
static_assert(std::endian::native == std::endian::little);
static_assert(std::numeric_limits<float>::is_iec559 && std::numeric_limits<float>::digits == 24);
static_assert(std::is_standard_layout_v<g::RetainedTriple_0> && sizeof(g::RetainedTriple_0) == 20);
static_assert(offsetof(g::RetainedTriple_0, valid_0) == 16);
static_assert(sizeof(StructuredBuffer<uint32_t>) == 16);
static_assert(sizeof(RWStructuredBuffer<uint32_t>) == 16);
static_assert(offsetof(StructuredBuffer<uint32_t>, data) == 0);
static_assert(offsetof(StructuredBuffer<uint32_t>, count) == 8);
static_assert(offsetof(RWStructuredBuffer<uint32_t>, data) == 0);
static_assert(offsetof(RWStructuredBuffer<uint32_t>, count) == 8);
static_assert(sizeof(g::GlobalParams_0) == 32 && offsetof(g::GlobalParams_0, inputs_0) == 0);
static_assert(offsetof(g::GlobalParams_0, outputs_0) == 16);
static_assert(std::is_standard_layout_v<g::KernelContext_0> && sizeof(g::KernelContext_0) == CONTEXT_SIZE);
OFFSETS
static_assert(std::is_same_v<decltype(&g::sirius_SLUGComputeMain_Group), void(*)(ComputeVaryingInput*, void*, void*)>);
static_assert(std::is_same_v<decltype(&g::sirius_SLUGComputeMain_Thread), void(*)(ComputeThreadVaryingInput*, void*, void*)>);
static_assert(sizeof(sirius_cpu_bridge::SLUG_program) / sizeof(uint32_t) == WORDS);
extern "C" bool sirius_native_SLUG_try_packet(
    const void* input, uint32_t input_bytes, void* output, uint32_t output_bytes,
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y,
    uint32_t grid_z, uint32_t groups) noexcept {
    if (!sirius_cpu_bridge::admit_packet(input, input_bytes, output, output_bytes,
            first_x, y, z, grid_x, grid_y, grid_z, groups, INPUT_ROW, OUTPUT_ROW,
            sirius_cpu_bridge::SLUG_program, WORDS))
        return false;
    g::GlobalParams_0 params{};
    params.inputs_0.data = const_cast<uint32_t*>(static_cast<const uint32_t*>(input));
    params.inputs_0.count = input_bytes / sizeof(uint32_t);
    params.outputs_0.data = static_cast<uint32_t*>(output);
    params.outputs_0.count = output_bytes / sizeof(uint32_t);
    for (uint32_t i = 0; i < groups; ++i) {
        ComputeVaryingInput varying{};
        varying.startGroupID.x = first_x + i;
        varying.startGroupID.y = y;
        varying.startGroupID.z = z;
        g::sirius_SLUGComputeMain_Group(&varying, nullptr, &params);
    }
    return true;
}
'''.replace('NAMESPACE', namespace).replace('CONTEXT_SIZE', str(context_size))
       .replace('OFFSETS', assertions).replace('INPUT_ROW', str(input_row))
       .replace('OUTPUT_ROW', str(output_row)).replace('WORDS', str(words)).replace('SLUG', slug))
    output = B / (slug + '.o')
    argv = [compiler, '-std=c++2c', '-O2', '-fPIC', '-DNDEBUG',
            '-Wall', '-Wextra', '-Wpedantic', '-Werror',
            '-isystem', '/opt/slang/include', '-isystem', str(B), '-I', str(W),
            '-march=x86-64', '-mtune=znver4', '-mfpmath=sse', '-fexcess-precision=standard',
            '-fno-fast-math', '-fno-associative-math', '-fno-unsafe-math-optimizations',
            '-frounding-math', '-ffp-contract=off', '-c', str(unit), '-o', str(output)]
    inputs.extend((original, reflection, header, canonical, unit))
    prepared.append({'stage': stage, 'argv': argv, 'output': str(output.relative_to(R)),
                     'complete_program_words': words, 'context_bytes': context_size})
record = {
    'objects': prepared,
    'inputs': {str(path.relative_to(R)): pin(path) for path in inputs},
    'target_fixes_per_stage': 'Exactly4 precise tokens removed under strictFP; exactly1 per-thread context zero-initialized. Generated namespace+4 uniquely named entries avoid ODR/C-linkage collisions. No other generated operations/control edited.',
    'target': 'Baseline x86-64 ISA/SSE2, tuned for znver4; no znver4 ISA requirement. Little-endian/IEEE32/pointer/layout assertions.',
    'scope': 'Prepared source only; no compiler, native execution, provider integration, FP equivalence or performance claim.',
}
(W / 'prepared-batch-inputs.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'stages': len(prepared), 'sealed_sources': len(record['inputs']),
                  'compiler_executed': False}))
