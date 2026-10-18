"""Read first control8 or all8 original native controls, without execution."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
from native_evidence import Context, native_case, timestamp_result, factored_routes, document, local, identity, pin, publish


def main():
    mode, bindings_arg, admission_arg, output_arg = sys.argv[1:]
    assert mode in ('timestamp', 'all8')
    context = Context(Path(bindings_arg).resolve(strict=True), Path(admission_arg).resolve(strict=True))
    assert context.value['phase'] == ('timestamp' if mode == 'timestamp' else 'all8')
    labels = ['control-8'] if mode == 'timestamp' else ['control-' + str(i) for i in range(1, 9)]
    assert set(context.value['cases']) == set(labels)
    cases = {label:native_case(context, label) for label in labels}
    assert len({tuple(case['windows_birth']) for case in cases.values()}) == len(labels)
    timestamp = timestamp_result(cases['control-8']['props'], 'candidate', context)
    rows = [{'label': label, 'case': cases[label]['owner']['expected_cases'][0],
             'accepted': cases[label]['accepted'], 'actual_windows_birth': cases[label]['windows_birth'],
             'wall_seconds': cases[label]['owner']['elapsed_seconds'],
             'peak_sampled_child_WorkingSet_bytes': cases[label]['owner']['peak_sampled_rss_bytes']}
            for label in labels]
    result = {'schema': 'sirius-native-packed-original-controls-acceptance-v1', 'passed': True,
              'source_revision': context.value['producer_source_revision'], 'control8_passed': True,
              'original_controls_passed': len(labels), 'controls': rows, 'timestamp': timestamp,
              'scope': 'Original native integrated controls with unchanged scientific/refusal/recovery assertions, same-child mapped SDK/AMD and terminal ownership. No full qualification, cold-frame, release or production adoption.'}
    if mode == 'all8':
        benefit_row = context.value['benefit_acceptance']
        benefit_path = local(benefit_row['path'])
        assert identity(benefit_path) == pin(benefit_row) == pin(context.files[benefit_row['path']])
        benefit = document(benefit_path)
        assert benefit['schema'] == 'sirius-native-packed-original22-benefit-v1'
        assert benefit['passed'] is True and benefit['audit_passed'] is True and benefit['benefit_passed'] is True
        assert benefit['source_revision'] == context.value['producer_source_revision']
        assert type(benefit['checks_total']) is int and benefit['checks_total'] == 22
        assert type(benefit['checks_passed']) is int and benefit['checks_passed'] == 22
        assert timestamp['readbacks'] == benefit['records']['a1']['readbacks']
        result['factored_endpoint_routes'] = factored_routes(cases['control-1'], context)
        result['all8_original_controls_and_both_complete_readers_passed'] = True
        result['benefit_acceptance'] = benefit_row
        result['readback_scope'] = '32 timestamp buffers exact across current native ABBA arms;28 native Factored buffers complete/hashed with unchanged in-test reference and refusal assertions. No fabricated equality to historical portable-provider bytes.'
    publish(Path(output_arg), result, context)


if __name__ == '__main__':
    main()
