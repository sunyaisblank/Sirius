"""Read the original affected control without executing or weakening its assertions."""
from pathlib import Path
import hashlib,json,xml.etree.ElementTree as ET
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'control-5'
D=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
action=D(W/'actions.json')['control-5'];owner=D(F/'owner.json')
assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0 and owner['stop_reason'] is None
assert owner['source']=={'revision':action['revision'],'status':''}
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
assert owner['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
assert owner['wall_seconds']<=300 and owner['peak_sampled_rss_kib']<=4194304 and owner['minimum_sampled_available_kib']>=2097152
before=D(F/'inputs-before.json');assert before==D(F/'inputs-after.json')
for name,v in before.items():assert pin(R/name)=={k:v[k] for k in ('bytes','sha256')},name
assert owner['declared_output_environment']=={'GTEST_OUTPUT':'xml:'+str(F/'gtest.xml'),'VK_DRIVER_FILES':action['runtime_icd'],'SIRIUS_VULKAN_DEVICE':'0','MESA_SHADER_CACHE_DISABLE':'true'}
assert owner['provider_mappings'] and all(row['owned_session']==owner['child_pid'] and row['maps'] and all(action['runtime_provider'] in line for line in row['maps']) for row in owner['provider_mappings'])
xml=ET.parse(F/'gtest.xml').getroot();assert int(xml.get('tests','-1'))==1
for node in [xml,*xml.findall('.//testsuite')]:
 for key in ('failures','errors','skipped','disabled'):assert int(node.get(key,'0'))==0
assert not any(xml.findall('.//'+tag) for tag in ('failure','error','skipped'))
cases=xml.findall('.//testcase');assert len(cases)==1
case=cases[0];assert case.get('classname')+'.'+case.get('name')==action['case']
assert case.get('status')=='run' and case.get('result')=='completed'
result={'passed':True,'case':action['case'],'owner':pin(F/'owner.json'),'XML':pin(F/'gtest.xml'),
 'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],
 'scope':action['scope'],'reader':pin(Path(__file__).resolve())}
with (W/'result.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'original_control5_passed':True,'wall_seconds':owner['wall_seconds']}))
