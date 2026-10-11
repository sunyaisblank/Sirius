"""Join the eight original native controls before permitting matched timing commands."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import importlib.util,json,re
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('native_comparison',WORK/'compare_matched.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
review=m.document(WORK/'controls-reader-v2-review.json')
assert review['accepted'] and m.whole_seal(Path(__file__))==review['reader_seal']
plan=m.document(WORK/'execution-plan.json')
head=plan['source_revision']
names=re.findall(r'^- `([A-Za-z]+Test\.[A-Za-z]+)`$',(WORK/'execution-preregistration.md').read_text(),re.M)
cases=m.document(WORK/'native-cases.json')['cases']
assert list(cases.values())==names and len(names)==8
records={};input_seals={};owner_births=set();test_births=set();reference_provider=None
for case,filter in cases.items():
 input_seals.update(m.validate_execution(case))
 folder=WORK/case/'run'
 owner=m.owner(folder,case);props=m.original_properties(folder,case);provider=m.provider(folder)
 if reference_provider is None:reference_provider=provider
 assert provider==reference_provider
 owner_birth=(owner['owner_pid'],owner['owner_start_time_utc'])
 test_birth=(owner['test_pid'],owner['test_start_time_utc'])
 assert owner_birth not in owner_births and test_birth not in test_births
 owner_births.add(owner_birth);test_births.add(test_birth)
 records[case]={'filter':filter,'normal_completed':True,'owner_birth':owner_birth,'test_birth':test_birth,
                'elapsed_seconds':owner['elapsed_seconds'],'peak_sampled_rss_bytes':owner['peak_sampled_rss_bytes'],
                'bridge_seal':m.whole_seal(WORK/case/'native-control-bridge.json'),
                'retirement_seal':m.whole_seal(WORK/case/'retirement-run.json')}
f=m.original_properties(WORK/'factored/run','factored')
t=m.original_properties(WORK/'timestamps/run','timestamps')
fr=m.readbacks(f,'factored');tr=m.readbacks(t,'timestamps');lr=m.legacy_readbacks(t)
ft=m.factored_timings(WORK/'factored/run');ts=m.timestamps(t)
arrays=m.document(WORK/'product-inspection.json')['arrays']
stage_arrays={'camera':'Camera','transport':'Transport','endpoint':'Endpoint','dense':'Dense',
              'initialize':'Initialize','ray_camera':'RayCamera','dopri_phase':'DopriPhase'}
bool_values={'false':0,'true':1}
assert t['fma_fp32_enabled'] in bool_values and t['preserves_fp32_signed_zero_inf_nan'] in bool_values
fma_enabled=bool_values[t['fma_fp32_enabled']]
signed_zero=bool_values[t['preserves_fp32_signed_zero_inf_nan']]
header=(WORK/'candidate/retained_kernels.h').read_text()
selected=m.selected_timestamp_modules(t,ft[0]['selected_stage_sha256'],arrays)
for route in ft:
 assert route['fma_enabled']==fma_enabled
 for stage,array in stage_arrays.items():
  suffix='Fp64Shader' if route['wide'] else 'Shader'
  if route['wide'] and stage in ('transport','endpoint'):
   availability=re.findall(r'inline constexpr bool k'+array+r'FmaAvailable = (true|false);',header)
   assert len(availability)==1
   if route['fma_enabled'] and signed_zero and availability[0]=='true':
    suffix='FmaShader'
  assert route['selected_stage_sha256'][stage]==arrays['k'+array+suffix]['sha256'], (stage,suffix)
for path in (Path(__file__),WORK/'compare_matched.py',WORK/'analysis-review.json',WORK/'controls-reader-v2-review.json',WORK/'execution-plan.json'):
 input_seals[str(path.relative_to(ROOT))]=m.whole_seal(path)
for name,expected in input_seals.items():assert m.whole_seal(ROOT/name)==expected,name
report={'accepted':True,'kind':'Eight original integrated native contract controls','source_revision':head,
        'reader':m.identity(Path(__file__)),
        'execution_plan_seal':m.whole_seal(WORK/'execution-plan.json'),'input_seals':input_seals,
        'all_eight_original_cases_normal_completed_zero_failure_error_skip_disabled':True,
        'records':records,'actual_provider_modules':reference_provider,
        'complete_factored_readbacks':fr,'complete_coupled_readbacks':tr,'complete_legacy_readbacks':lr,
        'timestamp_selected_modules':selected,'ordered_available_markers':51,
        'factored_routes':ft,'modern_coupled':ts,'source_product_and_bridge_seals_joined':True,
        'matched_timing_or_adoption_claimed':False,'full_scientific_frame_release_qualification_claimed':False,
        'scope':'Eight original numerical/reference/refusal/recovery/governor and ordered-stream controls, with exact source/provider/product/bridge/retirement joins. Complete same-provider buffer identity and frozen44-condition benefit require separate matched acceptance.'}
with (WORK/'controls-analysis.json').open('x') as stream:stream.write(json.dumps(report,indent=2)+'\n')
print(json.dumps({'accepted':True,'original_controls':8,'factored_readbacks':len(fr),'coupled_readbacks':len(tr),'ordered_markers':51,'matched_adoption_pending':True}))
