"""Read the one frozen comparison; complete-call benefit decides adoption."""
import json, math, statistics
from native_contracts import WORK, document, read_case, sha
bindings=document(WORK/'execution-bindings.json')
freeze=document(WORK/'retention-gate-frozen.json')
assert freeze['comparator_sha256']==sha(WORK/'compare_native.py')
assert freeze['issue_readback_sha256']==sha(WORK/'issue-readback.json')
sequence=document(WORK/'native-matched-sequence.json')
assert sequence['completed'] and sequence['commands']==[{'label':label,'case':'timestamps'} for label in ('a1','b1','b2','a2')]
records={label:read_case('matched',label,'timestamps',bindings) for label in ('a1','b1','b2','a2')}
for label,record in records.items():
    terminal=document(WORK/('terminal-matched-'+label+'-timestamps-accepted.json'));assert terminal['pass'] and terminal['record']['artifacts']==record['artifacts']
    original=records['a1']
    assert record['timestamps']['readbacks']==original['timestamps']['readbacks']
    assert record['timestamps']['identity']==original['timestamps']['identity'] and record['modules']==original['modules']
    expected='0' if label.startswith('a') else '1'
    assert all(m['submission_fence']==expected for family in ('legacy','modern') for r in record['timestamps'][family]['repeats'] for m in r['markers'])
checks=[]
def compare(direction,family,scope,a,b,threshold):
    ma=statistics.median(a);mb=statistics.median(b);assert math.isfinite(ma) and ma>0 and math.isfinite(mb) and mb>0
    improvement=100*(ma-mb)/ma
    checks.append({'direction':direction,'family':family,'scope':scope,'baseline_samples_ms':a,'candidate_samples_ms':b,'baseline_median_ms':ma,'candidate_median_ms':mb,'improvement_percent':improvement,'minimum_percent':threshold,'pass':improvement>=threshold})
for direction,a_label,b_label in [('forward','a1','b1'),('reverse','a2','b2')]:
    for family in ('legacy','modern'):
        a=records[a_label]['timestamps'][family]['repeats'];b=records[b_label]['timestamps'][family]['repeats']
        compare(direction,family,'complete_call_ms',[r['complete_call_ms'] for r in a],[r['complete_call_ms'] for r in b],5)
        stages={m['stage'] for r in a for m in r['markers']}
        assert stages=={m['stage'] for r in b for m in r['markers']}
        for stage in sorted(stages):
            compare(direction,family,stage+' inclusive Dispatch',[math.fsum(m['total_ms'] for m in r['markers'] if m['stage']==stage) for r in a],[math.fsum(m['total_ms'] for m in r['markers'] if m['stage']==stage) for r in b],-5)
assert len(checks)==20 and sum(q['scope']=='complete_call_ms' for q in checks)==4
result={'source_revision':bindings['candidate_revision'],'baseline_revision':bindings['baseline_revision'],'pass':all(q['pass'] for q in checks),'checks':checks,'failed_checks':[q for q in checks if not q['pass']],'readback_identity_module_and_completion_prerequisites_pass':True,'scope':'Exactlyone original native A1/B1/B2/A2 experiment, sourcebound defaultbinary32 Schwarzschild12-row modern cohort and legacyKerr/flat. No frame/fullscience/platform/release claim.'}
output=WORK/'native-comparison.json';assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'pass':result['pass'],'checks':len(checks),'failed':len(result['failed_checks']),'complete_calls':[q for q in checks if q['scope']=='complete_call_ms']},indent=2))
