"""Challenge the default timestamp reader with consequential evidence defects."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import copy, json, hashlib
from native_contracts import ROOT, WORK, document, sha, xml_properties, timestamp_result
path=WORK/'native/baseline/baseline/timestamps/gtest.xml'
props=xml_properties(path,'timestamps');positive=timestamp_result(props,'baseline')
assert len(positive['readbacks'])==32 and len(positive['selected_modules'])==7
cases={
 'logical_fma_disabled':lambda p:p.__setitem__('fma_fp32_enabled','0'),
 'signed_zero_control_missing':lambda p:p.__setitem__('preserves_fp32_signed_zero_inf_nan','false'),
 'wrong_product_mode':lambda p:p.__setitem__('product_mode','fp64'),
 'wrong_driver':lambda p:p.__setitem__('selected_driver','invented'),
 'wrong_timestamp_period':lambda p:p.__setitem__('timestamp_period_ns','1'),
 'missing_complete_buffer':lambda p:p.pop('coupled_readback_10'),
 'wrong_complete_buffer_size':lambda p:p.__setitem__('coupled_readback_10',p['coupled_readback_10'].replace('bytes=','bytes=9')),
 'missing_selected_shader':lambda p:p.pop('selected_module_6'),
 'wrong_selected_transport_payload':lambda p:p.__setitem__('selected_module_1',p['selected_module_1'].replace('sha256=','sha256=0')),
 'missing_ordered_marker':lambda p:p.pop('device_observation_50'),
 'unavailable_ordered_marker':lambda p:p.__setitem__('device_observation_0',p['device_observation_0'].replace('available0=1','available0=0')),
 'wrong_ordered_marker_phase':lambda p:p.__setitem__('device_observation_0',p['device_observation_0'].replace('phase=legacy_repeat0','phase=legacy_repeat1')),
 'nonfinite_whole_call':lambda p:p.__setitem__('legacy_repeat0_complete_call_ms','NaN'),
}
results=[]
for name,mutate in cases.items():
 bad=copy.deepcopy(props);mutate(bad);assert bad!=props
 try:timestamp_result(bad,'baseline')
 except (AssertionError,ValueError,KeyError) as e:results.append({'case':name,'rejected':True,'exception_type':type(e).__name__})
 else:raise AssertionError('Reader accepted invalid evidence: '+name)
result={'pass':True,'positive_original_native_baseline':True,'positive_xml':{'path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,'sha256':sha(path)},'negative_controls':results,'comparator_sha256':sha(WORK/'compare_native.py'),'scope':'Finite source/evidence reader controls only; no additional scientific dispatch or timing result and no predicate/threshold change.'}
(WORK/'native-reader-controls.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'pass':True,'positive':1,'negative_rejections':len(results)}))
