"""Decide the one preregistered complete-call ABBA trial."""
import json, math, statistics
from software_contracts import WORK, document, timestamp_result, xml_properties, sha
freeze=document(WORK/'retention-gate-frozen.json')
assert freeze['reader_sha256']==sha(WORK/'software_contracts.py')
assert freeze['comparator_sha256']==sha(WORK/'compare_software.py')
records={}
for label in ('a1','b1','b2','a2'):
    producer='baseline' if label.startswith('a') else 'candidate'
    folder=WORK/label; owner=document(folder/'owner.json')
    assert owner['passed'] and owner['status']=='terminal' and owner['child_exit']==0
    assert not owner['stop_reason'] and not owner['cleanup_errors'] and owner['provider_mappings']
    assert owner['remaining_owned_processes']==[] and all(b['absent'] for b in owner['observed_births'])
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['source']=={'revision':freeze['candidate_revision'],'status':''}
    assert owner['argv']==[freeze['executables'][producer]['path'],'--gtest_filter=RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults']
    props=xml_properties(folder/'gtest.xml','timestamps')
    records[label]=timestamp_result(props,producer)
    original=records['a1']
    assert records[label]['readbacks']==original['readbacks']
    assert records[label]['identity']==original['identity']
    old=original['selected_modules']; new=records[label]['selected_modules']
    assert {stage for stage in old if new[stage]!=old[stage]}==(set() if producer=='baseline' else {'transport','endpoint'})
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
result={'candidate_revision':freeze['candidate_revision'],'baseline_revision':freeze['baseline_revision'],
    'pass':all(check['pass'] for check in checks),'checks':checks,'failed_checks':[q for q in checks if not q['pass']],
    'readback_identity_module_and_completion_prerequisites_pass':True,
    'scope':'One fixed software default A1/B1/B2/A2 original complete observer. Four complete-call >=5% and sixteen stage >=-5% gates; no pooled orders/retries, frame, science or qualification claim.'}
path=WORK/'software-comparison.json'; assert not path.exists()
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'pass':result['pass'],'checks':20,'failed':len(result['failed_checks']),
    'complete_calls':[q for q in checks if q['scope']=='complete_call_ms']},indent=2))
