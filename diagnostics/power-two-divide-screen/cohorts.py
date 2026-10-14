"""Classify source-proven divisor groups; never estimate dynamic frequency."""
from pathlib import Path
import collections
import importlib.util
import json
import sys
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('issue99_cohorts', 'src/sirius/kernels/retained_program.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
results = []
for builder, target in [('build_transport_program', 2557), ('build_schwarzschild_transport_program', 427)]:
    program = getattr(module, builder)(parallel=True, with_plan=True)
    plan = program['source_plan']
    ops, writes = plan['operations'], plan['writes']
    stats, groups = collections.Counter(), []
    for layer, (start, end) in enumerate(zip(program['layer_offsets'], program['layer_offsets'][1:])):
        for base in range(start, end, 8):
            members = []
            for index in range(base, min(base + 8, end)):
                node = writes[index]
                op = ops[node]
                if op[0] != 5:
                    continue
                denominator = ops[op[2]]
                kind = 'unknown'
                if denominator[0] == 8:
                    kind = 'constructed_power_of_two'
                elif denominator[0] == 0:
                    word = denominator[1] & 0x7fffffff
                    if word and word < 0x7f800000 and (
                            word & 0x7fffff == 0 or word < 0x800000 and word & (word - 1) == 0):
                        kind = 'constant_unit' if word == 0x3f800000 else 'constant_power_of_two'
                members.append({'write_index': index, 'source_node': node,
                    'denominator_source_node': op[2], 'denominator_source': denominator, 'kind': kind})
            if not members:
                continue
            all_known = all(x['kind'] != 'unknown' for x in members)
            stats['divide_simd_groups'] += 1
            stats['all_source_known_power_or_unit_groups'] += all_known
            stats['some_source_known_power_or_unit_groups'] += any(x['kind'] != 'unknown' for x in members)
            if all_known or any(x['write_index'] == target for x in members):
                groups.append({'layer': layer, 'write_indices': [base, min(base + 8, end)],
                    'all_division_lanes_power_or_unit': all_known, 'members': members})
    results.append({'builder': builder, 'stats': dict(stats), 'interesting_groups': groups})
output = {'source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
          'scope': 'Source dependency and SIMD8 cohort classification, using current normal Transport stride16 and actual LLVM vector width8. No dynamic frequency, successful arithmetic count, exclusive cost or speed claim. Earlier exact-unit/validity guards can remove calls before PB32Divide.',
          'corrected_initial_diagnostic_filter': 'Use compiled write index, not original source-DAG node id, for 2557/427 loci.',
          'results': results}
Path(__file__).with_name('source-cohorts.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps([{'builder': r['builder'], **r['stats']} for r in results]))
