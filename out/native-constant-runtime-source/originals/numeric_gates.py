"""Exact original 22 numerical comparisons only; no provenance/runtime admission.

Caller must first independently admit four genuine native records with the
original internal three repeats, full raw readbacks, markers, modules and owners.
This function does not supply or assert that prerequisite proof.
"""
import math
import statistics


def evaluate_validated_records(records, freeze):
    checks=[]
    def compare(direction,family,scope,a,b,threshold):
        ma=statistics.median(a); mb=statistics.median(b)
        assert math.isfinite(ma) and ma>0 and math.isfinite(mb) and mb>0
        improvement=100*(ma-mb)/ma
        checks.append({'direction':direction,'family':family,'scope':scope,'baseline_samples_ms':a,
            'candidate_samples_ms':b,'baseline_median_ms':ma,'candidate_median_ms':mb,
            'improvement_percent':improvement,'minimum_percent':threshold,'pass':improvement>=threshold})
    for direction,a_label,b_label in [('forward','a1','b1'),('reverse','a2','b2')]:
        for family in ('legacy','modern'):
            a=records[a_label][family]['repeats']; b=records[b_label][family]['repeats']
            compare(direction,family,'complete_call_ms',[r['complete_call_ms'] for r in a],[r['complete_call_ms'] for r in b],5)
            stages={m['stage'] for r in a for m in r['markers']}
            assert stages=={m['stage'] for r in b for m in r['markers']}
            for stage in sorted(stages):
                compare(direction,family,stage+' inclusive Dispatch',
                    [math.fsum(m['total_ms'] for m in r['markers'] if m['stage']==stage) for r in a],
                    [math.fsum(m['total_ms'] for m in r['markers'] if m['stage']==stage) for r in b],-5)
    assert len(checks)==20 and sum(q['scope']=='complete_call_ms' for q in checks)==4
    assert freeze['observed_instance_api_minimum_improvement_percent'] == 0
    accounting_checks = []
    for direction, a_label, b_label in [('forward','a1','b1'),('reverse','a2','b2')]:
        a = records[a_label]['instance_api_duration_ms']; b = records[b_label]['instance_api_duration_ms']
        assert math.isfinite(a) and a > 0 and math.isfinite(b) and b > 0
        improvement = 100 * (a - b) / a
        accounting_checks.append({'direction':direction,'scope':'observed-instance Create and all11 original numerical API windows',
            'baseline_ms':a,'candidate_ms':b,'improvement_percent':improvement,'minimum_percent':0,'pass':improvement>=0})
    return {
        'all22_numeric_pass': all(q['pass'] for q in checks + accounting_checks),
        'original20_pass': all(q['pass'] for q in checks),
        'provisioning_accounting_pass': all(q['pass'] for q in accounting_checks),
        'checks': checks,
        'accounting_checks': accounting_checks,
        'scope': 'Pure original numerical predicates only; no native evidence or adoption acceptance.',
    }
