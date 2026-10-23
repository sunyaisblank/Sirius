"""Compile the generated shader with two audited target fixes; never execute it."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-slang-bridge/native-object'
original = R / 'bin/linux-gcc/native-cpu-slang-bridge/primitive-lowering/dopri_phase.cpp'
generated = original.read_text()
assert generated.count('precise float ') == 4
assert generated.count('KernelContext_0 kernelContext_10;') == 1
fixed = generated.replace('precise float ', 'float ')
fixed = fixed.replace('KernelContext_0 kernelContext_10;',
                      'KernelContext_0 kernelContext_10{};')
B.mkdir(parents=True, exist_ok=False)
header = B / 'dopri_phase.native.hpp'
header.write_text(fixed)
unit = B / 'typed_dopri.cpp'
unit.write_text('''#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <dopri_phase.native.hpp>
static_assert(sizeof(void*) == 8 && sizeof(float) == 4);
static_assert(std::endian::native == std::endian::little);
static_assert(std::numeric_limits<float>::is_iec559 && std::numeric_limits<float>::digits == 24);
static_assert(std::is_standard_layout_v<RetainedTriple_0> && sizeof(RetainedTriple_0) == 20);
static_assert(offsetof(RetainedTriple_0, valid_0) == 16);
static_assert(sizeof(StructuredBuffer<uint32_t>) == 16);
static_assert(sizeof(RWStructuredBuffer<uint32_t>) == 16);
static_assert(offsetof(StructuredBuffer<uint32_t>, data) == 0);
static_assert(offsetof(StructuredBuffer<uint32_t>, count) == 8);
static_assert(offsetof(RWStructuredBuffer<uint32_t>, data) == 0);
static_assert(offsetof(RWStructuredBuffer<uint32_t>, count) == 8);
static_assert(sizeof(GlobalParams_0) == 32 && offsetof(GlobalParams_0, inputs_0) == 0);
static_assert(offsetof(GlobalParams_0, outputs_0) == 16);
static_assert(std::is_standard_layout_v<KernelContext_0> && sizeof(KernelContext_0) == 12592);
static_assert(offsetof(KernelContext_0, retainedStatus_0) == 8);
static_assert(offsetof(KernelContext_0, retainedScratch_0) == 12);
static_assert(std::is_same_v<decltype(&ComputeMain_Group), void(*)(ComputeVaryingInput*, void*, void*)>);
static_assert(std::is_same_v<decltype(&ComputeMain_Thread), void(*)(ComputeThreadVaryingInput*, void*, void*)>);
''')
compiler = str(Path(shutil.which('g++-14')).resolve(strict=True))
output = B / 'dopri_phase.o'
args = [compiler, '-std=c++2c', '-O2', '-fPIC', '-DNDEBUG',
        '-Wall', '-Wextra', '-Wpedantic', '-Werror',
        '-isystem', '/opt/slang/include', '-isystem', str(B),
        '-march=znver4', '-mfpmath=sse', '-fexcess-precision=standard',
        '-fno-fast-math', '-fno-associative-math', '-fno-unsafe-math-optimizations',
        '-frounding-math', '-ffp-contract=off', '-c', str(unit), '-o', str(output)]
result = subprocess.run(args)
def pin(p):
    data = p.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
record = {'argv': args, 'exit': result.returncode, 'original_CPP': pin(original),
          'fixed_generated_header': pin(header), 'typed_assertions_unit': pin(unit),
          'scope': 'One C++26 compile-only ABI/object screen. Generated code and installed prelude are isolated as compiler-produced system headers; owned typed ABI assertions retain strict warnings. FP policy is explicit and applies to all code. No object execution, provider integration, runtime/FP equivalence or benefit claim.',
          'target_fixes': 'Exactly four shader-only precise tokens removed under enforced no-fast-math/no-reassociation/FP-contract-off, and one per-invocation KernelContext aggregate zero initialization. No other generated operation/control change.'}
if output.exists():
    record['actual_object'] = {'path': str(output.relative_to(R)), **pin(output)}
(W / 'compile-result.json').write_text(json.dumps(record, indent=2) + '\n')
raise SystemExit(result.returncode)
