"""Evaluate the unchanged preregistered 44-condition complete-call gate."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,importlib.util,json,subprocess
from datetime import datetime,timezone
w=Path(__file__).resolve().parent;r=w.parent.parent
def doc(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=r,text=True,timeout=20).strip()
spec=importlib.util.spec_from_file_location('native_comparison',w/'compare_matched.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
plan=doc(w/'execution-plan.json');head=plan['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
for name,expected in plan['input_seals'].items():assert m.whole_seal(r/name)==expected,name
review=doc(w/'analysis-review.json');assert review['accepted']
for name,expected in review['analysis_seals'].items():assert m.whole_seal(w/name)==expected,name
c=doc(w/'controls-analysis.json');d=doc(w/'matched-native-comparison.json')
assert c['accepted'] and c['source_revision']==d['live_source_revision']==head
assert c['execution_plan_seal']==m.whole_seal(w/'execution-plan.json')
assert c['all_eight_original_cases_normal_completed_zero_failure_error_skip_disabled'] and c['source_product_and_bridge_seals_joined']
assert c['actual_provider_modules']==d['selected_provider_modules']
assert all(d[k] is True for k in ['all_28_factored_complete_readbacks_equal','all_11_coupled_complete_readbacks_equal','all_51_original_marker_order_and_availability_equal','all_eight_original_owners_completed','bridge_pre_post_source_and_product_seals_checked_here','all_21_legacy_complete_readbacks_equal','all_seven_timestamp_selected_modules_joined_to_actual_producer'])
assert c['complete_legacy_readbacks']==d['legacy_readbacks'] and c['complete_factored_readbacks']==d['factored_readbacks'] and c['complete_coupled_readbacks']==d['coupled_readbacks']
for name,expected in c['input_seals'].items():assert m.whole_seal(r/name)==expected,name
for record in d['input_files'].values():
 assert m.whole_seal(r/record['path'])=={k:record[k] for k in ('bytes','sha256')}
assert d['observer']==m.identity(w/'compare_matched.py')
bridge_records=[]
actions=doc(w/'native-cases.json')['actions'];assert len(actions)==16
for action in actions:
 m.validate_execution(action)
 p=w/action/'native-control-bridge.json';b=doc(p)
 bridge_records.append({'path':str(p.relative_to(r)),'sha256':sha(p),'producer':b['producer'],'action':action})
ci=doc(w/'ci-completed.json')
assert ci['databaseId']==38159862169 and ci['headSha']==head and ci['status']=='completed' and ci['conclusion']=='success'
assert all(any(j['name']==name and j['conclusion']=='success' for j in ci['jobs']) for name in ('integration-no-render','integration-macos-no-render','integration-windows-no-render'))
checks = []
def check(direction, scope, metric, value, bound, relation):
    pct = value['candidate_percent_change']
    passed = pct < bound if relation == '<' else pct <= bound
    checks.append({'direction':direction,'scope':scope,'metric':metric,**value,'required_percent_change':f'{relation} {bound}', 'passed':passed})
for direction, q in d['directions'].items():
    for metric in ['warm13_wait_ms','warm13_total_ms']:
        check(direction,'default Endpoint',metric,q['factored_routes'][0][metric],-5,'<=')
        check(direction,'excluded wide Endpoint',metric,q['factored_routes'][1][metric],5,'<=')
    for metric in ['device_ms','host_completion_ms']:
        check(direction,'complete coupled sum',metric,q['coupled_marker_sums'][metric],-5,'<=')
    check(direction,'complete coupled sum','complete_call_ms',q['coupled_complete_calls']['complete_call_ms'],-5,'<=')
    for index, repeat in enumerate(q['coupled_per_repeat_calls']):
        check(direction,f'coupled repeat {index}','complete_call_ms',repeat['complete_call_ms'],0,'<')
    for repeat, metrics in q['coupled_per_repeat_markers'].items():
        for metric in ['device_ms','host_completion_ms']:
            check(direction,repeat,metric,metrics[metric],0,'<=')
    for stage in ['transport','dopri_phase']:
        for metric in ['device_ms','host_completion_ms','total_ms']:
            check(direction,('targeted Transport' if stage=='transport' else 'excluded '+stage),metric,q['coupled_stage_sums'][stage][metric],5,'<=')
failed=[x for x in checks if not x['passed']]
assert len(checks)==44
report={'kind':'Root joined original native bound-selector retention decision',
 'utc':datetime.now(timezone.utc).isoformat(),'candidate_revision':head,'baseline_revision':plan['baseline_revision'],
 'current_source_clean_committed_pushed':True,'current_complete_input_seals_exact':True,
 'current_complete_input_seal_count':len(plan['input_seals']),'all_eight_original_controls_passed':True,
 'all_eight_matched_owners_completed':True,'all_28_factored_and_11_coupled_buffers_equal':True,
 'all_21_legacy_buffers_equal':True,'all_seven_timestamp_modules_joined':True,'all_51_markers_ordered_available':True,'same_loaded_provider_and_device':True,
 'input_reports':{name:sha(w/name) for name in ('controls-analysis.json','matched-native-comparison.json','execution-preregistration.md','execution-preregistration-readback.json','execution-plan.json','ci-completed.json','analysis-review.json')},
 'bridge_records':bridge_records,'retention_checks':checks,'failed_conditions':failed,
 'submit_and_wait_drift':{key:value['coupled_stage_sums'] for key,value in d['directions'].items()},
 'decision':'reject benefit; preserve trial and restore accepted baseline in pushed Git' if failed else 'finite original integrated numerical and benefit gate passed; broader rendering/qualification remains open',
 'full_frame_rerun_authorized_by_this_result':False,'full_scientific_release_qualification_claimed':False,
 'limits':'Finite original controls and same-provider complete-buffer identities, not exhaustive numerical proof. Stage spans include synchronization/scheduling, and cold cache is unverified. Existing <=5% Transport drift conditions remain required although Transport is also targeted. Raw ignored evidence remains local only.'}
with (w/'root-retention-decision.json').open('x') as stream:stream.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'decision':report['decision'],'conditions':len(checks),'failed_conditions':len(failed),'failed_scopes':sorted(set(x['scope'] for x in failed))}))
