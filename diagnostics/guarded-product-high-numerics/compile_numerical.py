"""Compile private arithmetic probes with original recipes and narrow FMul review."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024**2, 64 * 1024**2))
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/guarded-product-high-numerics'
family, producer, mode = sys.argv[1:]
assert producer in ['stock', 'candidate']
assert family in ['residual', 'limb']
assert mode in (['core', 'rte32'] if family == 'residual' else
                ['integer', 'normal', 'normal_transport', 'normal_endpoint'])
assert family != 'residual' or producer == 'candidate' or mode == 'core'
include = stage / producer / 'kernels'
source = work / ('residual_probe.slang' if family == 'residual' else 'limb_arithmetic_probe.slang')
stem = stage / f'{family}-{producer}-{mode}'
raw, optimized, assembly, output = [stem.with_suffix(suffix)
    for suffix in ['.raw.spv', '.optimized.spv', '.spvasm', '.spv']]
assert not any(path.exists() for path in [raw, optimized, assembly, output])
cache = (root / 'bin/linux-gcc/CMakeCache.txt').read_text()
def tool(name):
    return re.search(r'^' + name + r':FILEPATH=(.*)$', cache, re.M)[1]
compiler, assembler, disassembler, validator, optimizer = [tool(name) for name in
    ['SIRIUS_SLANGC', 'SIRIUS_SPIRV_AS', 'SIRIUS_SPIRV_DIS', 'SIRIUS_SPIRV_VAL', 'SIRIUS_SPIRV_OPT']]
definitions = ['-DSIRIUS_RETAINED_PORTABLE=1']
normal = family == 'limb' and mode != 'integer'
if normal:
    definitions.append('-DSIRIUS_RETAINED_NORMAL_SUM32=1')
if mode == 'normal_transport':
    definitions.append('-DSIRIUS_RETAINED_PARALLEL_TRANSPORT=1')
if mode == 'normal_endpoint':
    definitions.append('-DSIRIUS_RETAINED_PARALLEL_ENDPOINT=1')
commands = []
def run(argv):
    commands.append(argv)
    subprocess.run(argv, check=True)
run([compiler, str(source), *definitions, '-I', str(include), '-O0', '-target', 'spirv',
     '-profile', 'spirv_1_5', '-entry', 'ComputeMain', '-stage', 'compute', '-o', str(raw)])
run([optimizer, '--target-env=vulkan1.2', '--ssa-rewrite', '--eliminate-dead-code-aggressive',
     '--preserve-bindings', '--preserve-interface', str(raw), '-o', str(optimized)])
run([disassembler, str(optimized), '-o', str(assembly)])
text = assembly.read_text()
spec = importlib.util.spec_from_file_location('original_retained_controls', root / 'scripts/build-retained-kernels.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)
multiplies = re.findall(r'(%\S+) = OpFMul (%\S+) (%\S+) (%\S+)', text)
assert re.findall(r'OpCapability (\S+)', text) == ['Shader']
assert set(re.findall(r'OpTypeInt (\d+) [01]', text)) == {'32'}
assert not re.search(r'\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision|OpFmaKHR)\b', text)
assert not re.search(r'\bOp(?:ConvertFTo[SU]|Convert[SU]ToF|FConvert|Dot|QuantizeToF16|VectorTimes\w+|MatrixTimes\w+)\b', text)
private_allowance = False
if producer == 'candidate':
    assert len(multiplies) == 1
    value, float_type, _, _ = multiplies[0]
    assert re.findall(r'OpTypeFloat (\d+)', text) == ['32']
    assert f'OpDecorate {value} NoContraction' in text
    bodies = [body for body in re.findall(r'%\S+ = OpFunction .*?OpFunctionEnd', text, re.S)
              if '= OpFMul ' in body]
    assert len(bodies) == 1 and re.match(r'%PB32ProductResidual(?:_\d+)? = OpFunction ', bodies[0])
    assert 'PB32HighModel' not in text
    if normal:
        try:
            builder.portable_normal_sum_controls(text)
        except ValueError as error:
            assert str(error) == 'portable normal sum changed its ordered three-operation body'
        else:
            raise AssertionError('original checker unexpectedly admitted extra FMul')
        filtered = re.sub(r'(?m)^.* = OpFMul .*\n', '', text)
        builder.portable_normal_sum_controls(filtered)
        private_allowance = True
    else:
        assert re.findall(r'= (OpF(?!unction)\S+) ', text) == ['OpFMul']
        private_allowance = True
elif normal:
    assert not multiplies
    text = builder.portable_normal_sum_controls(text)
else:
    assert 'OpTypeFloat' not in text and not multiplies
requires_rte32 = normal or mode == 'rte32'
if producer == 'candidate' and requires_rte32:
    entry = re.search(r'OpEntryPoint GLCompute (%\S+)', text)[1]
    text = text.replace('OpCapability Shader', 'OpCapability Shader\n               OpCapability RoundingModeRTE', 1)
    local = re.search(r'^.*OpExecutionMode ' + re.escape(entry) + r' LocalSize.*$', text, re.M)[0]
    text = text.replace(local, local + '\n               OpExecutionMode ' + entry + ' RoundingModeRTE 32', 1)
assembly.write_text(text)
run([assembler, '--target-env', 'spv1.5', str(assembly), '-o', str(output)])
run([validator, '--target-env', 'vulkan1.2', str(output)])
final = stem.with_suffix('.final.spvasm')
run([disassembler, str(output), '-o', str(final)])
actual = final.read_text()
assert re.findall(r'OpCapability (\S+)', actual) == (['Shader', 'RoundingModeRTE'] if requires_rte32 else ['Shader'])
assert re.findall(r'OpExecutionMode %\S+ LocalSize (\d+ \d+ \d+)', actual) == ['64 1 1']
assert bool(re.search(r'OpExecutionMode %\S+ RoundingModeRTE 32', actual)) == requires_rte32
assert not re.search(r'\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision|OpFmaKHR)\b', actual)
assert len(re.findall(r'= OpFMul ', actual)) == (1 if producer == 'candidate' else 0)
record = lambda path: {'path': str(path.relative_to(root)), 'bytes': path.stat().st_size,
                       'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
report = {'family': family, 'producer': producer, 'mode': mode, 'commands': commands,
    'actual_module': record(output), 'actual_final_disassembly': record(final),
    'include_directory': str(include.relative_to(root)), 'candidate_header': record(include / 'portable_binary32.h'),
    'complete_actual_module_validated': True, 'normal_body_checker_unchanged': True,
    'private_guarded_fmul_allowance': private_allowance,
    'filtered_text_assembled_or_executed': False, 'required_rte32': requires_rte32,
    'scope': 'Actual complete private arithmetic module; added guarded FMul is not an adopted production policy.'}
(work / f'compiled-{family}-{producer}-{mode}.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: report[key] for key in ['family', 'producer', 'mode', 'actual_module', 'required_rte32']}))
