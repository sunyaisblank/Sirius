"""Apply the issue83 finite gate exactly once to the frozen matched sequence."""
from pathlib import Path
import json, math
from native_contracts import WORK, document, read_case

bindings = document(WORK / 'execution-bindings.json')
controls = document(WORK / 'native-controls-analysis.json')
assert controls['pass'] and len(controls['results']) == 9
sequence = document(WORK / 'native-matched-sequence.json')
assert sequence['completed'] and len(sequence['results']) == 4
records = {label: read_case('matched', label, 'timestamps', bindings) for label in ('a1', 'b1', 'b2', 'a2')}
reference = records['a1']['timestamps']
for record in [*records.values(), controls['results'][-1]]:
    assert record['timestamps']['readbacks'] == reference['readbacks']
    assert record['timestamps']['identity'] == reference['identity']
    assert record['modules'] == records['a1']['modules']
checks = []
def check(name, baseline, candidate, rule, threshold):
    assert rule in ('improvement_at_least', 'strict_improvement', 'no_regression', 'regression_at_most')
    assert math.isfinite(baseline) and math.isfinite(candidate) and baseline > 0
    improvement = 100 * (1 - candidate / baseline)
    passed = improvement >= threshold if rule == 'improvement_at_least' else candidate < baseline if rule == 'strict_improvement' else candidate <= baseline if rule == 'no_regression' else improvement >= -threshold
    checks.append({'name': name, 'baseline_ms': baseline, 'candidate_ms': candidate,
                   'improvement_percent': improvement, 'rule': rule, 'threshold_percent': threshold, 'pass': passed})
for direction, a, b in [('forward', 'a1', 'b1'), ('reverse', 'a2', 'b2')]:
    before, after = records[a]['timestamps'], records[b]['timestamps']
    for metric in ('device_ms', 'host_completion_ms', 'complete_call_ms'):
        check(direction + '.legacy.aggregate.' + metric, before['legacy']['sums'][metric], after['legacy']['sums'][metric], 'improvement_at_least', 5)
        check(direction + '.modern.aggregate.' + metric, before['modern']['sums'][metric], after['modern']['sums'][metric], 'regression_at_most', 5)
    for repeat in range(3):
        for metric, rule in [('complete_call_ms', 'strict_improvement'), ('device_ms', 'no_regression'), ('host_completion_ms', 'no_regression')]:
            check(f'{direction}.legacy.repeat{repeat}.{metric}', before['legacy']['repeats'][repeat][metric], after['legacy']['repeats'][repeat][metric], rule, 0)
    # The final combined Endpoint+Dense timestamp cannot be split into independent
    # stages. Keep that original non-Transport family combined for drift checks.
    for stage in ('endpoint', 'endpoint+dense'):
        for metric in ('device_ms', 'host_completion_ms'):
            def family(value):
                return math.fsum(m[metric] for r in value['legacy']['repeats'] for m in r['markers'] if m['stage'] == stage)
            check(direction + '.legacy.nontransport.' + stage + '.' + metric, family(before), family(after), 'regression_at_most', 5)
    for stage in ('endpoint', 'dopri_phase'):
        for metric in ('device_ms', 'host_completion_ms'):
            def family(value):
                return math.fsum(m[metric] for r in value['modern']['repeats'] for m in r['markers'] if m['stage'].split(':', 1)[0] == stage)
            check(direction + '.modern.nontransport.' + stage + '.' + metric, family(before), family(after), 'regression_at_most', 5)
result = {'candidate_revision': bindings['candidate_revision'], 'baseline_revision': bindings['baseline_revision'],
          'records': records, 'checks': checks, 'failed_checks': [v for v in checks if not v['pass']],
          'numerical_stream_provider_owner_prerequisites_pass': True,
          'retention_gate_pass': all(v['pass'] for v in checks),
          'scope': 'One finite original native A1/B1/B2/A2 comparison. Engineering 5% gate, no confidence, full-frame, scientific-estate or release claim.'}
assert len(checks) == 46
output = WORK / 'native-retention-comparison.json'
assert not output.exists()
output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'checks': len(checks), 'failed_checks': len(result['failed_checks']), 'retention_gate_pass': result['retention_gate_pass'], 'failed_names': [v['name'] for v in result['failed_checks']]}))
