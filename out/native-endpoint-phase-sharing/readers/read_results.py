"""Preserve original C8/22 predicates; bind both genuinely different providers."""
from pathlib import Path
from types import SimpleNamespace
import hashlib,importlib.util,json,sys
R=Path.cwd();W=Path(__file__).resolve().parent;RUN=R/'out/native-endpoint-phase-sharing/runtime';VALID=R/'out/native-endpoint-phase-sharing/provider-history'
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
old=load('original_timestamp_authority',W/'native_evidence.py');numeric=load('original22',W/'numeric_gates.py');old.WORK=W
def doc(p):return old.document(p)
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
actions=doc(RUN/'actions.json');providers=doc(W.parent/'providers.json')['providers']
native=dict(zip(old.STAGES,['k'+p+'Shader' for p in old.PREFIXES]));native['endpoint']='kEndpointNativeOwnedShader';native['dopri_phase']='kDopriPhaseNativeOwnedShader'
portable=dict(zip(old.STAGES,['kCameraPortableShader','kTransportPortableNormalSumShader','kEndpointPortableNormalSumShader','kDensePortableShader','kInitializePortableShader','kRayCameraPortableShader','kDopriPhasePortableNormalSumShader']))
routes={'baseline':portable,'candidate':native};old.selected_names=lambda context,arm,wide=False:routes[arm]
files={'baseline':['camera_portable','transport_portable_normal_sum','endpoint_portable_normal_sum','dense_portable','initialize_portable','ray_camera_portable','dopri_phase_portable_normal_sum'],'candidate':['camera','transport','endpoint','dense','initialize','ray_camera','dopri_phase']}
arrays={};hardware={}
for arm,label in [('baseline','baseline-poison-v21'),('candidate','candidate-default')]:
    for stage,stem in zip(old.STAGES,files[arm]):
        p=(R/doc(W.parent/'compiler-result.json')['modules'][{'endpoint':'Endpoint','dopri_phase':'DopriPhase'}[stem]]['module'] if arm=='candidate' and stem in ['endpoint','dopri_phase'] else R/'bin/linux-gcc/src/sirius/backend/retained'/('retained_'+stem+'.spv'))
        arrays[routes[arm][stage]]=pin(p)
    lines=(VALID/label/'stdout.log').read_text().splitlines();info,=[json.loads(x) for x in lines if x.startswith('{') and 'resident_driver_path' in x]
    assert info['kind']=='software' and info['driver']=='llvmpipe' and info['driver_id']==13 and info['vendor_id']==65541 and info['device_id']==0
    assert info['denorm32']==int(arm=='candidate') and info['RTE32']==1 and 'LLVM 21.1.8' in info['device']
    expected_provider=providers[arm];assert info['resident_driver_path']==expected_provider['historical_library'] and pin(Path(expected_provider['library']))==expected_provider['library_identity']
    props=doc(VALID/label/'float-controls-properties.json');assert all(x['values']['shaderDenormPreserveFloat32']==int(arm=='candidate') and x['values']['shaderRoundingModeRTEFloat32']==1 for x in props)
    creation,=doc(VALID/label/'device-features.json')['device_creations'];fma=any(x['fma16_32_64'][1]==1 for x in creation['feature_chain'])
    hardware[arm]={'selected_device':info['device'],'selected_driver':info['driver']+': '+info['driver_info'],'retained_scalar_route':'native binary32' if arm=='candidate' else 'portable with guarded normal sums','timestamp_valid_bits':'64','timestamp_period_ns':'1','fma_fp32_enabled':str(int(fma)),'preserves_fp32_signed_zero_inf_nan':str(props[0]['values']['shaderSignedZeroInfNanPreserveFloat32'])}
assert hardware['baseline']['selected_device']==hardware['candidate']['selected_device']
def admit(label):
    folder=RUN/label;owner=doc(folder/'owner.json');action=actions[label];arm=action['arm']
    assert owner['passed'] and owner['status']=='terminal' and owner['child_exit']==0 and owner['stop_reason'] is None
    assert owner['cleanup_errors']==owner['remaining_owned_processes']==[]
    assert owner['source']=={'revision':action['revision'],'status':''} and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    assert owner['wall_seconds']<=300 and owner['peak_sampled_rss_kib']<=4194304 and owner['minimum_sampled_available_kib']>=2097152
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
    before=doc(folder/'inputs-before.json');assert before==doc(folder/'inputs-after.json')
    for key,value in before.items():
        p=Path(key) if key.startswith('/') else R/key;assert pin(p)=={k:value[k] for k in ['bytes','sha256']},key
    env=owner['declared_output_environment'];assert env['VK_DRIVER_FILES']==providers[arm]['ICD'] and env['SIRIUS_VULKAN_DEVICE']=='0' and env['MESA_SHADER_CACHE_DISABLE']=='true' and env['GTEST_OUTPUT']=='xml:'+str(folder/'gtest.xml')
    assert owner['provider_mappings'] and all(x['owned_session']==owner['child_pid'] and x['maps'] and all(providers[arm]['library'] in line for line in x['maps']) for x in owner['provider_mappings'])
    props=old.xml_properties(folder/'gtest.xml',label)
    ctx=SimpleNamespace(value={'native_identity':hardware[arm]},arms={x:{'arrays':arrays} for x in ['baseline','candidate']})
    return owner,old.timestamp_result(props,arm,ctx)
control_owner,control=admit('control-8');mode,=sys.argv[1:]
if mode=='control':
    result={'passed':True,'record':control,'owner':pin(RUN/'control-8/owner.json'),'scope':'Original C8,51 markers54 executed-module joins32 full readbacks and unchanged assertions, under two ownership modules on the actual native route enabled by implemented support. No benefit/adoption/frame/release acceptance.'};output=W/'control-8-acceptance.json'
elif mode=='benefit':
    assert doc(W/'control-8-acceptance.json')['passed'];records={};owners={};last=control_owner['child_start_ticks'];births=set()
    for label in ['a1','b1','b2','a2']:
        owner,record=admit(label);assert owner['boot_id']==control_owner['boot_id'];birth=(owner['child_pid'],owner['child_start_ticks']);assert birth not in births and birth[1]>last;births.add(birth);last=birth[1]
        assert record['readbacks']==control['readbacks']
        # Provider implementation now honestly supports native32. Only its
        # labelled driver version and derived implementation route may differ.
        assert {k:v for k,v in record['identity'].items() if k not in ['selected_driver','retained_scalar_route']}=={k:v for k,v in control['identity'].items() if k not in ['selected_driver','retained_scalar_route']}
        expected=set(old.STAGES) if actions[label]['arm']=='baseline' else set();assert {s for s in old.STAGES if record['selected_modules'][s]!=control['selected_modules'][s]}==expected
        records[label]=record;owners[label]={'owner':pin(RUN/label/'owner.json'),'birth':birth,'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib']}
    decision=numeric.evaluate_validated_records(records,{'observed_instance_api_minimum_improvement_percent':0});checks=decision['checks']+decision['accounting_checks'];assert len(checks)==22
    result={'passed':decision['all22_numeric_pass'],'audit_passed':True,'benefit_passed':decision['all22_numeric_pass'],'checks_total':22,'checks_passed':sum(x['pass'] for x in checks),'fixed_order':['a1','b1','b2','a2'],'internal_repeats_per_slot':3,'extra_retries_or_warmups':0,'records':records,'owners':owners,**decision,'scope':'One original fixed ABBA22 on the same ELF/physicalCPU/LLVM21 with implemented provider precision support. Original4 whole-call>=5%,16 inclusive stages>=-5%,2 provisioning-inclusive API>=0%;32 cross-arm raw readback digests exact. No pooling/subtraction/cache/governor/input/count/threshold changes or retries. Remaining7 and cold frame required after22 PASS.'};output=W/'benefit.json'
else:raise ValueError(mode)
assert not output.exists();output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['passed','audit_passed','checks_total','checks_passed','scope'] if k in result}))
