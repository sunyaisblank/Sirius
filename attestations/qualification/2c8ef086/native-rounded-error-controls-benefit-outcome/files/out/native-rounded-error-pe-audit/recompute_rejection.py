import hashlib
import json
import re
import xml.etree.ElementTree as ET
from decimal import Decimal, getcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
getcontext().prec = 45
path = ROOT / 'out/native-rounded-error-ABBA-evidence/native-benefit.json'
report = json.loads(path.read_text())
assert report['audit_passed'] and not report['benefit_passed'] and report['checks_passed'] == 20
props, pins = {}, {}
for label in ['a1', 'b1', 'b2', 'a2']:
    xml = ROOT / 'out/native-rounded-error-runtime-abba/native' / label / 'gtest.xml'
    raw = xml.read_bytes()
    pairs = [(item.get('name'), item.get('value')) for item in ET.fromstring(raw).findall('.//property')]
    assert len(dict(pairs)) == len(pairs)
    props[label] = dict(pairs)
    pins[str(xml.relative_to(ROOT))] = {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def median(values):
    assert len(values) == 3 and all(value > 0 and value.is_finite() for value in values)
    return sorted(values)[1]


def calls(values, prefix):
    return [Decimal(values[f'{prefix}{i}_complete_call_ms']) for i in range(3)]


def stage_sums(values, prefix):
    result = {}
    for index in range(51):
        fields = dict(pair.split('=', 1) for pair in values[f'device_observation_{index}'].split(';'))
        phase = fields['phase']
        if not phase.startswith(prefix):
            continue
        repeat = int(phase.removeprefix(prefix))
        assert 0 <= repeat < 3
        result.setdefault(fields['stage'], [Decimal(0)] * 3)[repeat] += Decimal(fields['total_ms'])
    return result


checks = []
for direction, baseline, candidate in [('forward', 'a1', 'b1'), ('reverse', 'a2', 'b2')]:
    for family, prefix, count in [('legacy', 'legacy_repeat', 3), ('modern', 'coupled_repeat', 5)]:
        pairs = [('complete_call_ms', calls(props[baseline], prefix), calls(props[candidate], prefix), 5)]
        left, right = stage_sums(props[baseline], prefix), stage_sums(props[candidate], prefix)
        assert set(left) == set(right) and len(left) == count
        pairs += [(stage + ' inclusive Dispatch', left[stage], right[stage], -5) for stage in sorted(left)]
        for scope, left_values, right_values, minimum in pairs:
            gain = Decimal(100) * (median(left_values) - median(right_values)) / median(left_values)
            original, = [item for item in report['checks'] if (item['direction'], item['family'], item['scope']) == (direction, family, scope)]
            passed = gain >= minimum
            assert passed == original['pass'] and abs(gain - Decimal(str(original['improvement_percent']))) < Decimal('0.00000001')
            checks.append({'direction': direction, 'family': family, 'scope': scope,
                           'raw_XML_improvement_percent': str(gain), 'minimum_percent': minimum, 'passed': passed})
    totals = []
    for label in [baseline, candidate]:
        values = props[label]
        assert Decimal(values['owner_provisioning_ms']) == 0
        totals.append(sum([Decimal(values[key]) for key in ['observed_create_ms', 'disabled_attempt_ms', 'disabled_sample_ms']]
                          + calls(values, 'legacy_repeat') + calls(values, 'coupled_repeat')))
    gain = Decimal(100) * (totals[0] - totals[1]) / totals[0]
    original, = [item for item in report['accounting_checks'] if item['direction'] == direction]
    assert (gain >= 0) == original['pass'] and abs(gain - Decimal(str(original['improvement_percent']))) < Decimal('0.00000001')
    checks.append({'direction': direction, 'scope': original['scope'], 'raw_XML_improvement_percent': str(gain),
                   'minimum_percent': 0, 'passed': gain >= 0})
assert len(checks) == 22 and sum(item['passed'] for item in checks) == 20
readbacks = None
for label, values in props.items():
    current = {key: value for key, value in values.items()
               if re.fullmatch(r'legacy_repeat\d+_readback_\d+|coupled_readback_\d+', key)}
    assert len(current) == 32 and (readbacks is None or current == readbacks)
    readbacks = current
seed = {item.get('name'): item.get('value') for item in ET.parse(ROOT / 'out/native-rounded-error-runtime-timestamp/native/control-8/gtest.xml').findall('.//property')}
assert all(seed[key] == value for key, value in readbacks.items())
assert not (ROOT / 'out/native-rounded-error-runtime-remaining').exists()
assert not (ROOT / 'out/native-rounded-error-all8-evidence').exists()
review = {'schema': 'sirius-native-rounded-error-root-raw-rejection-check-v1', 'passed': True,
          'producer_source_revision': report['source_revision'], 'raw_XML': pins,
          'implementation': 'Standalone Decimal45 XML property arithmetic; no original reader or numeric helper imported.',
          'all22_outcomes_and_improvements_reproduced': True, 'checks': checks,
          'same_complete32_readback_properties_across_control8_and_four_slots': True,
          'remaining7_not_materialized': True, 'benefit_passed': False,
          'scope': 'Separate ROOT arithmetic/raw-data check; original complete951-input ownership/provider/readback audit remains the mature reader evidence. No new workload, retry, sample, adoption or full qualification.'}
out = WORK / 'root-raw-rejection-check.json'
assert not out.exists()
out.write_text(json.dumps(review, indent=2) + '\n')
decision = {'schema': 'sirius-native-rounded-error-root-disposition-v1', 'decision': 'REJECT',
            'producer_source_revision': report['source_revision'], 'canonical_source_revision': '2c8ef08693d46f57d5dea14efdbddd4874dac684',
            'audit_passed': True, 'control8_passed': True, 'benefit_passed': False, 'checks_passed': 20, 'checks_total': 22,
            'whole_calls_passed': 2, 'whole_calls_total': 4, 'stages_passed': 16, 'stages_total': 16,
            'inclusive_API_passed': 2, 'inclusive_API_total': 2,
            'failed_conditions': [item for item in checks if not item['passed']],
            'remaining7_state': 'not materialized or run', 'source_adopted': False,
            'unchanged_retry_or_pooling_or_resampling': False, 'frame_run_performed': False,
            'cause_of_forward_reverse_difference': 'Unproved; no noise/cache/thermal/causal attribution.',
            'scope': 'Exact error shortcut is not admitted for adoption after sole originalABBA benefit refusal. Preserve negative outcome; pursue a distinct supported cost premise.'}
(ROOT / 'out/native-rounded-error-ci-publication/root-rejection.json').write_text(json.dumps(decision, indent=2) + '\n')
print(json.dumps({'raw_recomputation_passed': True, 'audit_passed': True, 'benefit_passed': False, 'checks_passed': 20, 'remaining7_not_run': True}))
