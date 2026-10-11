"""Read completed original scientific controls without launching work."""
from pathlib import Path
import hashlib, json, sys, xml.etree.ElementTree as ET
WORK=Path(__file__).resolve().parent
CASES={'candidate-rk':'RetainedComputeTest.JointRkStagesRetainCriticalIncrementsAndEmbeddedError',
       'candidate-mixed':'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal',
       'candidate-schwarzschild':'RetainedComputeTest.SchwarzschildStagesPreserveIndependentFieldsAndGeneralFallback',
'candidate-cache-primer':'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
'candidate-rk-cached':'RetainedComputeTest.JointRkStagesRetainCriticalIncrementsAndEmbeddedError',
'candidate-mixed-cached':'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal',
'candidate-schwarzschild-cached':'RetainedComputeTest.SchwarzschildStagesPreserveIndependentFieldsAndGeneralFallback'}

def seal(p):
    r=p.read_bytes();return {'bytes':len(r),'sha256':hashlib.sha256(r).hexdigest()}

def main(name):
    assert name in CASES
    folder=WORK/name;owner=json.loads((folder/'owner.json').read_text())
    assert owner['passed'] and owner['status']=='terminal' and owner['child_exit']==0 and owner['stop_reason'] is None
    assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert all(r['absent'] for r in owner['observed_births'])
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['provider_mappings']
    before=json.loads((folder/'inputs-before.json').read_text());after=json.loads((folder/'inputs-after.json').read_text());assert before==after
    actions=json.loads((WORK/'action-versions'/(name+'.json')).read_text());a=actions[name]
    assert owner['argv']==a['argv'] and owner['source']=={'revision':a['revision'],'status':''} and owner['limits']==a['limits']
    for path,row in before.items():
        p=WORK/'action-versions'/(name+'.json') if path==str((WORK/'actions.json').relative_to(WORK.parents[1])) else Path(row['resolved_path'])
        assert seal(p)=={k:row[k] for k in ('bytes','sha256')},str(p)
    x=ET.fromstring((folder/'gtest.xml').read_bytes());assert x.tag=='testsuites' and x.attrib['tests']=='1'
    for node in [x,*x.findall('.//testsuite')]:
        assert all(int(node.attrib.get(k,'0'))==0 for k in ('failures','errors','skipped','disabled'))
    assert not any(x.findall('.//'+t) for t in ('failure','error','skipped'))
    cases=x.findall('.//testcase');assert len(cases)==1
    case=cases[0];assert case.attrib['classname']+'.'+case.attrib['name']==CASES[name] and case.attrib['status']=='run' and case.attrib['result']=='completed'
    extra = {}
    if name == 'candidate-cache-primer':
        import software_contracts as contracts
        contracts.ROOT = WORK.parents[1]; contracts.WORK = WORK
        observed = contracts.timestamp_result(contracts.xml_properties(folder/'gtest.xml', 'timestamps'), 'candidate')
        accepted = json.loads((WORK/'numerical-control-readback.json').read_text())['current_result']
        assert observed['readbacks'] == accepted['readbacks'] and len(observed['readbacks']) == 32
        assert observed['identity'] == accepted['identity'] and observed['selected_modules'] == accepted['selected_modules']
        assert owner['cache_before'] == {} and owner['cache_after']
        extra = {'complete_raw_digests_exact':32,'available_ordered_markers':51,'owned_cache_preparation_only':True,'cache_hit_unproved':True}
    result={'pass':True,**extra,'case':CASES[name],'xml':seal(folder/'gtest.xml'),'owner':seal(folder/'owner.json'),'action_version':seal(WORK/'action-versions'/(name+'.json')),'input_seals':len(before),'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'scope':'Original test body/assertions retained, one completedcase zeroF/E/skips/disabled and clean bounded owner/current input seals. Per-case LoadKernel/full-buffer properties only where original test actually publishes them; no invented raw comparisons, timing benefit or final qualification.'}
    (folder/'control-readback.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':
    assert len(sys.argv)==2;main(sys.argv[1])
