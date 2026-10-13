"""Compare bound ordinary complete-module code; do not infer runtime benefit."""
from pathlib import Path
import hashlib
import json
import re
import sys

sys.dont_write_bytecode = True
work = Path(__file__).resolve().parent
variants = ['endpoint-normal', 'endpoint-core', 'transport-normal',
            'dopri_phase-normal', 'dense-core']


def frame_stride(text, report):
    body = re.search(r'^define ptr @cs_co_variant\([^\n]*\{\n(.*?)^\}',
                     text, re.M | re.S).group(1)
    calls = re.findall(r'^\s*(%[\w.$-]+) = call ptr @coro_malloc\(i32 (%[\w.$-]+)\)$',
                       body, re.M)
    assert len(calls) == 1
    allocated, count = calls[0]
    definition = re.search(r'^\s*' + re.escape(count) + r' = (mul i32 (%[\w.$-]+), (\d+))$',
                           body, re.M)
    assert definition and definition.group(2) == '%12'
    stride = int(definition.group(3))
    offsets = re.findall(r'^\s*(%[\w.$-]+) = mul i32 %17, ' + str(stride) + r'$', body, re.M)
    assert len(offsets) == 1
    extended = re.findall(r'^\s*(%[\w.$-]+) = sext i32 ' + re.escape(offsets[0]) + r' to i64$',
                          body, re.M)
    assert len(extended) == 1
    assert re.search(r'getelementptr i8, ptr %[\w.$-]+, i64 ' + re.escape(extended[0]) + r'$', body, re.M)
    assert report['llvm_resume_frame_dereferenceable_bytes'] == [stride]
    return {'bytes': stride, 'allocation_size_ssa': count,
            'allocation_multiplication': definition.group(1),
            'allocation_call': '@coro_malloc(i32 ' + count + ')',
            'frame_index_offset_ssa': offsets[0],
            'frame_index_extended_ssa': extended[0],
            'allocation_stride_index_and_resume_bounds_agree': True}


def read(producer, variant):
    folder = work / f'capture-{producer}-{variant}'
    owner = json.loads((folder / 'owner.json').read_text())
    report = json.loads((folder / 'actual-code-review.json').read_text())
    assert owner['passed'] and owner['child_exit'] == 0 and not owner['stop_reason']
    export = report['ordinary_ir_export']
    raw = (folder / export['path']).read_bytes()
    assert len(raw) == export['bytes'] and hashlib.sha256(raw).hexdigest() == export['sha256']
    functions = {f['name']: f['bytes'] for f in report['jit_defined_functions']}
    return {'case': folder.name, 'actual_module': report['actual_loaded_module'],
            'optimized_ir': export, 'executable_section_bytes': report['jit_executable_section_bytes'],
            'resume_function_bytes': functions['cs_co_variant.resume'],
            'coroutine_frame': frame_stride(raw.decode(), report),
            'instrumented_capture_wall_seconds': owner['wall_seconds'],
            'instrumented_peak_sampled_rss_kib': owner['peak_sampled_rss_kib']}


if __name__ == '__main__':
    comparisons = []
    for variant in variants:
        baseline, candidate = [read(p, variant) for p in ['stock', 'candidate']]
        comparisons.append({'variant': variant, 'stock': baseline, 'candidate': candidate,
            'executable_section_growth_bytes': candidate['executable_section_bytes'] - baseline['executable_section_bytes'],
            'executable_section_growth_percent': 100 * (candidate['executable_section_bytes'] / baseline['executable_section_bytes'] - 1),
            'coroutine_frame_growth_bytes': candidate['coroutine_frame']['bytes'] - baseline['coroutine_frame']['bytes']})
    target = work / 'jit-comparison.json'
    assert not target.exists()
    result = {'source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
        'completed_single_dispatch_compiler_canaries': 10,
        'retained_arithmetic_exercised': False,
        'runtime_benefit_verdict': None,
        'production_adopted': False,
        'scope': 'Bound ordinary complete-module LLVM/JIT footprint comparison only. Instrumented capture duration is not retained runtime cost or a coldness claim.',
        'comparisons': comparisons}
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps([{k: row[k] for k in ['variant', 'executable_section_growth_bytes',
        'executable_section_growth_percent', 'coroutine_frame_growth_bytes']} for row in comparisons], indent=2))
