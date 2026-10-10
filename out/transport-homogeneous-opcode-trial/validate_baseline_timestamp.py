from pathlib import Path
import hashlib,json,math,subprocess,xml.etree.ElementTree as ET
from native_contracts import WORK,ROOT,xml_properties,CASES,MODERN,LEGACY,fields,number,sha,document,MARKER_MS
folder=WORK/'baseline-timestamp-controls/timestamps'
owner=document(folder/'owner.json');head=document(WORK/'test-first-source.json')['source_revision']
assert owner['source_revision']==head and owner['returncode']==0 and owner['stop_reason'] is None
assert owner['observer_error'] is None and not owner['cleanup_errors']
assert owner['owned_birth_absent'] and owner['owned_group_absent'] and owner['source_and_inputs_unchanged']
assert owner['timeout_seconds']==90 and owner['maximum_rss_mib']==4096
for path,record in owner['whole_input_seals'].items():
 p=Path(path);assert p.stat().st_size==record['bytes'] and sha(p)==record['sha256'],path
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
props=xml_properties(folder/'tests.xml','timestamps')
archive=ROOT/'attestations/native-vulkan/b9e8a1a/transport-critical-ready-rejected'
ref=archive/'controls/timestamps/tests.xml';manifest=document(archive/'manifest.json')
record=next(q for q in manifest['files'] if q['path']=='controls/timestamps/tests.xml')
assert ref.stat().st_size==record['bytes'] and sha(ref)==record['sha256']
reference={p.get('name'):p.get('value') for p in ET.parse(ref).findall('.//property')}
keys={f'legacy_repeat{r}_readback_{i}' for r in range(3) for i in range(7)}|{f'coupled_readback_{i}' for i in range(11)}
assert {k for k in props if '_readback_' in k}==keys
assert all(props[key]==reference[key] for key in keys)
assert {k for k in props if k.startswith('device_observation_')}=={f'device_observation_{i}' for i in range(51)}
order=[(f'legacy_repeat{r}',stage) for r in range(3) for stage in LEGACY]
order += [(f'coupled_repeat{r}',f'{stage}:rows={rows}') for r in range(3) for stage,rows in MODERN]
bits=int(props['timestamp_valid_bits']);period=number(props['timestamp_period_ns']);assert bits==64 and period==1
for i,pair in enumerate(order):
 marker=fields(props[f'device_observation_{i}']);assert set(marker)=={'phase','stage','begin','end','available0','available1',*MARKER_MS}
 assert (marker['phase'],marker['stage'])==pair and marker['available0']==marker['available1']=='1'
 begin,end=int(marker['begin']),int(marker['end']);assert 0<=begin<2**64 and 0<=end<2**64
 nums={k:number(marker[k]) for k in MARKER_MS}
 assert math.isclose(nums['device_ms'],((end-begin)&((1<<bits)-1))*period/1e6,rel_tol=1e-9,abs_tol=1e-9)
 assert math.isclose(nums['host_completion_ms'],nums['host_submit_ms']+nums['host_wait_ms'],rel_tol=1e-9,abs_tol=1e-9)
 assert math.isclose(nums['total_ms'],sum(nums[k] for k in ('pipeline_ms','setup_ms','host_completion_ms','cleanup_ms')),rel_tol=1e-9,abs_tol=1e-9)
result={'name':'timestamps','case':CASES['timestamps'],'source_revision':head,'returncode':0,'failures':0,'errors':0,'skips':0,'elapsed_seconds':owner['elapsed_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'properties':props,'owner_sha256':sha(folder/'owner.json'),'xml_sha256':sha(folder/'tests.xml'),'all32_same_provider_readbacks_exact':True,'all51_original_markers_complete':True,'postprocessor_sha256':sha(Path(__file__)),'reference':{'path':str(ref.relative_to(ROOT)),**record},'interpretation_correction':'Executed observer selected pre-observation304 XML and raised KeyError only after complete child/owner/input/XML success. No scientific execution repeated. This separate postprocessor joins retained b6 accepted same-provider complete observations to their existing archived manifest. Reference is postprocessing authority, not a retrospectively asserted prelaunch input.'}
assert not (folder/'result.json').exists();(folder/'result.json').write_text(json.dumps(result,indent=2)+'\n')
(WORK/'baseline-timestamp-controls-complete.json').write_text(json.dumps({'source_revision':head,'existing_tests':1,'new_test_registrations':0,'results':[result],'full_qualification':False},indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('source_revision','elapsed_seconds','all32_same_provider_readbacks_exact','all51_original_markers_complete')}))
