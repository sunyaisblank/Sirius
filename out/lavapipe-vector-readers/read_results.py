"""Admit Linux receipts, then reuse the unchanged original C8 and22 predicates."""
from pathlib import Path
from types import SimpleNamespace
import hashlib,importlib.util,json,os,sys
R=Path.cwd();W=Path(__file__).resolve().parent;RUN=R/'out/lavapipe-vector-runtime'
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
old=module('original_control_format',R/'out/lavapipe-vector-control-source/native_evidence.py')
numeric=module('original22_predicates',R/'out/lavapipe-vector-control-source/numeric_gates.py')
def doc(p):return old.document(p)
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def pin(p):return {'bytes':p.stat().st_size,'sha256':sha(p)}
actions=doc(RUN/'actions.json')
corpus=doc(RUN/'candidate-corpus64/corpus-acceptance.json')
device=corpus['provider'];assert device['kind']=='software' and device['driver']=='llvmpipe' and device['driver_id']==13
assert device['RTE32']==1 and device['denorm32']==0
old.WORK=W
names=dict(zip(old.STAGES,['kCameraPortableShader','kTransportPortableNormalSumShader','kEndpointPortableNormalSumShader','kDensePortableShader','kInitializePortableShader','kRayCameraPortableShader','kDopriPhasePortableNormalSumShader']))
old.selected_names=lambda context,arm,wide=False:dict(names)
files=['retained_camera_portable.spv','retained_transport_portable_normal_sum.spv','retained_endpoint_portable_normal_sum.spv','retained_dense_portable.spv','retained_initialize_portable.spv','retained_ray_camera_portable.spv','retained_dopri_phase_portable_normal_sum.spv']
arrays={names[stage]:pin(R/'bin/linux-gcc/src/sirius/backend/retained'/name) for stage,name in zip(old.STAGES,files)}

def admit(label):
 p=RUN/label;owner=doc(p/'owner.json');a=actions[label]
 assert owner['passed'] and owner['status']=='terminal' and owner['child_exit']==0 and owner['stop_reason'] is None
 assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
 assert owner['source']=={'revision':a['revision'],'status':''} and owner['source_unchanged'] and owner['whole_inputs_unchanged']
 assert owner['limits']==a['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
 assert owner['wall_seconds']<=300 and owner['peak_sampled_rss_kib']<=4194304 and owner['minimum_sampled_available_kib']>=2097152
 assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
 before=doc(p/'inputs-before.json');assert before==doc(p/'inputs-after.json')
 for key,value in before.items():
  path=Path(key) if key.startswith('/') else R/key
  assert pin(path)=={k:value[k] for k in ['bytes','sha256']},key
 env=owner['declared_output_environment']
 assert env['VK_DRIVER_FILES']==a['runtime_icd'] and env['SIRIUS_VULKAN_DEVICE']=='0' and env['MESA_SHADER_CACHE_DISABLE']=='true'
 assert env['GTEST_OUTPUT']=='xml:'+str(p/'gtest.xml')
 assert owner['provider_mappings']
 for item in owner['provider_mappings']:
  assert item['owned_session']==owner['child_pid'] and item['maps'] and all(a['runtime_provider'] in x for x in item['maps'])
 assert pin(Path(a['runtime_provider']))=={k:before[a['runtime_provider'].removeprefix(str(R)+'/')][k] for k in ['bytes','sha256']}
 props=old.xml_properties(p/'gtest.xml',label)
 return owner,props

control_owner,control_props=admit('control-8')
assert control_props['selected_device']==device['device']
assert control_props['selected_driver']==device['driver']+': '+device['driver_info']
assert control_props['retained_scalar_route']=='portable with guarded normal sums'
hardware={k:control_props[k] for k in ['selected_device','selected_driver','retained_scalar_route','timestamp_valid_bits','timestamp_period_ns','fma_fp32_enabled','preserves_fp32_signed_zero_inf_nan']}
assert hardware['fma_fp32_enabled'] in ['0','false']
ctx=SimpleNamespace(value={'native_identity':hardware},arms={a:{'arrays':arrays} for a in ['baseline','candidate']})
control=old.timestamp_result(control_props,'candidate',ctx)
mode,=sys.argv[1:]
if mode=='control':
 result={'schema':'sirius-private-lavapipe-original-C8-v1','passed':True,'record':control,'owner':pin(RUN/'control-8/owner.json'),
  'provider':pin(Path(actions['control-8']['runtime_provider'])),'scope':'Original C8 Linux software provider admission,51 markers54 executed-module joins32 complete readback digests, unchanged assertions. No ABBA/adoption/frame/release acceptance.'}
 output=W/'control-8-acceptance.json'
elif mode=='benefit':
 assert doc(W/'control-8-acceptance.json')['passed']
 records={};owners={};last=control_owner['child_start_ticks'];births=set()
 for label in ['a1','b1','b2','a2']:
  owner,props=admit(label);assert owner['boot_id']==control_owner['boot_id']
  birth=(owner['child_pid'],owner['child_start_ticks']);assert birth not in births and birth[1]>last
  births.add(birth);last=birth[1]
  arm='baseline' if label in ['a1','a2'] else 'candidate'
  record=old.timestamp_result(props,arm,ctx)
  assert record['readbacks']==control['readbacks'] and record['identity']==control['identity']
  assert record['selected_modules']==control['selected_modules']
  records[label]=record;owners[label]={'owner':pin(RUN/label/'owner.json'),'birth':birth,'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib']}
 decision=numeric.evaluate_validated_records(records,{'observed_instance_api_minimum_improvement_percent':0})
 checks=decision['checks']+decision['accounting_checks'];assert len(checks)==22
 result={'schema':'sirius-private-lavapipe-original22-v1','passed':decision['all22_numeric_pass'],'audit_passed':True,
  'benefit_passed':decision['all22_numeric_pass'],'checks_total':22,'checks_passed':sum(c['pass'] for c in checks),
  'fixed_order':['a1','b1','b2','a2'],'internal_repeats_per_slot':3,'extra_retries_or_warmups':0,'records':records,'owners':owners,
  **decision,'scope':'One original fixed matched software-provider trial,4 whole-call>=5%,16 inclusive stages>=-5%,2 provisioning-inclusive API>=0%. No pooling, subtraction, changes to inputs, caches, governor, thresholds or retries. Remaining7 controls and original cold frame required after22 PASS.'}
 output=W/'benefit.json'
else:raise ValueError(mode)
assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['passed','audit_passed','benefit_passed','checks_total','checks_passed','scope'] if k in result}))
