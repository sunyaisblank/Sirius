"""Exercise the real full-module declaration counterexample and retained policy."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import subprocess

root = Path.cwd()
output = root / 'out/portable-divide64-policy-controls.json'
assert not output.exists()
commit = 'da2338b6f2a803c9aa06801a3699163b25f80d6f'
source = ('attestations/software-vulkan/4bcd987/divide64-compiler-capability-order/'
          'files/bin/linux-gcc/portable-divide64-screen/modules/candidate-endpoint-normal.spvasm')
data = subprocess.check_output(['git', 'show', commit + ':' + source])
assert len(data) == 422604 and hashlib.sha256(data).hexdigest() == '2f7c83a973259e5becbacf572209e2160c7b3f7a13bc265ff765cdfbca59d753'
raw = data.decode()
builder = root / 'out/portable-divide64-screen-capability-order/candidate/scripts/build-retained-kernels.py'
before = hashlib.sha256(builder.read_bytes()).hexdigest()
spec = importlib.util.spec_from_file_location('full_module_policy', builder)
policy = importlib.util.module_from_spec(spec); spec.loader.exec_module(policy)
results = []
def check(name, text, expected, candidate=True):
    error = None
    try:
        policy.portable_normal_sum_controls(text, diagnostic_divide64=candidate)
        accepted = True
    except ValueError as exc:
        accepted = False; error = str(exc)
    assert accepted is expected, (name, expected, error)
    results.append({'name': name, 'expected_acceptance': expected, 'actual_acceptance': accepted, 'error': error,
                    'complete_fixture_sha256': hashlib.sha256(text.encode()).hexdigest()})
pattern = r'(?m)^([ \t]*)OpCapability (\S+)[ \t]*$'
assert re.findall(r'OpCapability (\S+)', raw) == ['Int64', 'Shader']
for order in (['Int64', 'Shader'], ['Shader', 'Int64']):
    names = iter(order)
    text = re.sub(pattern, lambda m: m[1] + 'OpCapability ' + next(names), raw)
    check('candidate-' + '-'.join(order), text, True)
    check('stock-rejects-' + '-'.join(order), text, False, candidate=False)
for cap in ('Int64', 'Int16'):
    check('extra-' + cap, raw.replace('OpCapability Shader', 'OpCapability Shader\n               OpCapability ' + cap, 1), False)
text, count = re.subn(r'OpTypeInt 64 0', 'OpTypeInt 16 0', raw, count=1)
assert count == 1; check('unsupported-integer-width', text, False)
check('unsupported-float-width', raw.replace('OpTypeFloat 32', 'OpTypeFloat 64', 1), False)
value = re.search(r'(%\S+) = OpF(?:Add|Sub) %\S+ ', raw)[1]
decoration = 'OpDecorate ' + value + ' NoContraction'
assert decoration in raw
check('relaxed-arithmetic', raw.replace(decoration, decoration + '\n               OpDecorate ' + value + ' RelaxedPrecision', 1), False)
text, count = re.subn(r'(?m)^[ \t]*' + re.escape(decoration) + r'[ \t]*\n', '', raw, count=1)
assert count == 1; check('missing-no-contraction', text, False)
text, count = re.subn(r' = OpFSub ', ' = OpFAdd ', raw, count=1)
assert count == 1; check('changed-normal-sum-order', text, False)
entry = re.search(r'OpEntryPoint GLCompute (%\S+)', raw)[1]
check('unqualified-float-mode', raw + '\nOpExecutionMode ' + entry + ' DenormPreserve 32\n', False)
assert len(results) == 12 and sum(r['actual_acceptance'] for r in results) == 2
assert hashlib.sha256(builder.read_bytes()).hexdigest() == before
output.write_text(json.dumps({'passed': True, 'controls': results, 'builder_sha256': before,
    'immutable_fixture_commit': commit, 'immutable_fixture_path': source,
    'immutable_fixture_sha256': hashlib.sha256(data).hexdigest(), 'full_original_module_retained': True,
    'compiler_or_device_executed': False, 'numerical_performance_or_qualification_claimed': False}, indent=2)+'\n')
print('PASS: two valid declaration orders and ten retained-policy negative controls')
