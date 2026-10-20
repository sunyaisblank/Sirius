"""Decide the sole original native A1/B1/B2/A2 trial from frozen receipts.

No process, provider, benchmark, resampling or retry is launched by this reader.
"""
import sys
sys.dont_write_bytecode = True
from datetime import datetime
from pathlib import Path
from native_evidence import Context, native_case, timestamp_result, publish, identity
from numeric_gates import evaluate_validated_records


def main():
    bindings_path, admission_path, output_path = map(Path, sys.argv[1:])
    context = Context(bindings_path.resolve(strict=True), admission_path.resolve(strict=True))
    assert set(context.value['cases']) == {'control-8', 'a1', 'b1', 'b2', 'a2'}
    assert context.value['phase'] == 'abba'
    initial = native_case(context, 'control-8')
    initial_data = timestamp_result(initial['props'], 'candidate', context)
    cases, records = {}, {}
    previous_end = None
    births = {tuple(initial['windows_birth'])}
    def instant(value):
        return datetime.fromisoformat(value.replace('Z', '+00:00'))
    for label in ('a1', 'b1', 'b2', 'a2'):
        case = native_case(context, label)
        assert tuple(case['windows_birth']) not in births
        births.add(tuple(case['windows_birth']))
        start, end = instant(case['owner']['test_start_time_utc']), instant(case['owner']['end_utc'])
        assert start < end
        if previous_end is not None: assert previous_end <= start
        else: assert instant(initial['owner']['end_utc']) <= start
        previous_end = end
        record = timestamp_result(case['props'], case['arm'], context)
        assert record['readbacks'] == initial_data['readbacks']
        assert record['identity'] == initial_data['identity']
        reference = records.get('a1' if case['arm'] == 'baseline' else 'b1')
        if reference is not None: assert record['selected_modules'] == reference['selected_modules']
        assert {k:v['array'] for k,v in record['selected_modules'].items()} == {k:v['array'] for k,v in initial_data['selected_modules'].items()}
        cases[label], records[label] = case, record
    decision = evaluate_validated_records(records, context.value)
    checks = decision['checks'] + decision['accounting_checks']
    assert len(checks) == 22
    passed = decision['all22_numeric_pass']
    result = {'schema': 'sirius-native-packed-original22-benefit-v1', 'passed': passed,
              'audit_passed': True, 'source_revision': context.value['producer_source_revision'],
              'benefit_passed': passed, 'checks_total': 22,
              'checks_passed': sum(q['pass'] for q in checks),
              'decision': 'ORIGINAL22_NATIVE_GATES_PASS' if passed else 'REJECT_NATIVE_PACKED_SCRATCH_COST',
              **decision, 'control8_passed': True,
              'fixed_order': ['a1', 'b1', 'b2', 'a2'], 'processes_per_slot': 1,
              'internal_repeats_per_family': 3, 'outer_repeats_or_warmup': 0,
              'records': records,
              'accepted_owners': {k:v['accepted'] for k,v in cases.items()},
              'actual_windows_births': {k:v['windows_birth'] for k,v in cases.items()},
              'source_readers': {name:identity(Path(__file__).parent/name) for name in ('native_evidence.py','numeric_gates.py','compare_native.py')},
              'prerequisites': {'same_native_identity': True, 'all32_raw_readbacks_equal': True,
                                'per_instance51markers54module_joins': True, 'actual_arrays_not_Linux_provider_maps': True},
              'scope': 'One fixed original native matched observer trial: 4 complete-call >=5%,16 inclusive stage >=-5%,2 complete API-window >=0%. No pooling, subtraction, warmup, retry, frame or scientific/release qualification. Remaining seven original controls and independent/root adoption decision are still required after benefit PASS.'}
    publish(output_path, result, context)


if __name__ == '__main__':
    main()
