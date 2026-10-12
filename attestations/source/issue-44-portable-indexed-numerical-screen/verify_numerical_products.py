"""Check four actual witness interfaces and retained indexed accumulator bodies."""
from pathlib import Path
import json
import re
import resource
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
from module_layout import module
root = Path.cwd()
work = Path(__file__).resolve().parent
build = root / 'bin/linux-gcc/portable-indexed-limb-accumulation'
old = root / 'attestations/software-vulkan/661e157/portable-limb-accumulation/files/bin/linux-gcc/portable-execution-attribution/limb-screen'
records = []
for mode in ('integer', 'normal', 'normal_transport', 'normal_endpoint'):
    before = module(old / ('arithmetic-baseline-' + mode + '.spv'))
    current = module(build / ('arithmetic-candidate-' + mode + '.spv'))
    for key in ('capabilities', 'integer_types', 'barriers', 'execution_modes',
                'descriptor_binding_layouts', 'no_contraction_operation_shapes'):
        assert before[key] == current[key], (mode, key)
    text = (build / ('arithmetic-candidate-' + mode + '.spvasm')).read_text()
    bodies = [m.group(0) for m in re.finditer(r'(%\S+) = OpFunction .*?OpFunctionEnd', text, re.S)
              if re.fullmatch(r'%RTAccumulateTerm(?:_\d+)?', m.group(1))]
    assert len(bodies) == 1
    body = bodies[0]
    calls = re.findall(r'= OpFunctionCall %\S+ (%\S+)', body)
    assert sum(v.startswith('%RetainedSum32') for v in calls) == 1
    assert sum(v.startswith('%RPUpAdd') for v in calls) == 1
    assert sum(v.startswith('%RPBounded') for v in calls) == 5
    assert body.count('OpLoopMerge') == 1 and 'DontUnroll' in body
    assert 'OpVariable %_ptr_Function__arr_uint_int_3 Function' in body
    assert re.search(r'OpULessThan %bool %\S+ %uint_3', body)
    assert len(current['descriptor_binding_layouts']) == 2
    records.append({'mode': mode, 'baseline': before, 'candidate': current,
                    'one_indexed_loop_one_sum_single_radius_five_bounds': True,
                    'accumulator_body': body})
result = {'passed': True, 'modes': records,
          'scope': 'Actual four portable witness products retain the existing interfaces/modes/controls and one indexed accumulator loop/body. Numerical, integrated source/stream, performance and qualification acceptance remain separate.'}
(work / 'numerical-products-check.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'modes': len(records)}))
