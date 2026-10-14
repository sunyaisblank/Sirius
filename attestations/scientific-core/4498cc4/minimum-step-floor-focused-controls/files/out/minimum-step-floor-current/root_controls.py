from pathlib import Path
import json,os,subprocess,sys,time,xml.etree.ElementTree as ET
work=Path(__file__).resolve().parent
controls=work/'controls'
actions=json.loads((controls/'actions.json').read_text())
commands=json.loads((controls/'commands.json').read_text())
assert json.loads((controls/'root-runtime-admission.json').read_text())['passed']
assert set(actions)==set(commands['order'])
assert all(not (controls/name).exists() for name in commands['order'])
environment=dict(os.environ,TMPDIR=str((work/'tmp').resolve()))
rows=[]
with (work/'root-controls-operator.log').open('xb') as log:
 for case in commands['order']:
  start=time.monotonic()
  result=subprocess.run(['/usr/bin/python3.12','-B',str(controls/'run_owned.py'),case],env=environment,stdout=log,stderr=subprocess.STDOUT)
  row={'case':case,'operator_exit':result.returncode,'operator_seconds':round(time.monotonic()-start,3),'accepted':False}
  try:
   folder=controls/case;owner=json.loads((folder/'owner.json').read_text());action=actions[case]
   assert result.returncode==0 and owner['passed'] and owner['status']=='terminal' and owner['child_exit']==0
   assert owner['source']=={'revision':'4498cc41168ad29fc4c111819a4cb7ccaa54a8bf','status':''}
   assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
   assert owner['limits']==action['limits'] and not owner['stop_reason'] and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
   assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
   assert json.loads((folder/'inputs-before.json').read_text())==json.loads((folder/'inputs-after.json').read_text())
   if action['software_provider']: assert owner['provider_mappings']
   xml=ET.parse(folder/('gtest.xml' if action['kind']=='test' else 'ctest.xml')).getroot()
   cases=xml.findall('.//testcase');assert len(cases)==1
   assert not any(xml.findall('.//'+tag) for tag in ['failure','error','skipped'])
   assert not any(int(node.get(key,'0')) for node in [xml,*xml.findall('.//testsuite')] for key in ['failures','errors','disabled','skipped'])
   testcase=cases[0]
   if action['kind']=='test':
    expected=action['root_required_xml']
    assert testcase.get('classname')==expected['suite'] and testcase.get('name')==expected['name']
    assert testcase.get('status')=='run' and testcase.get('result')=='completed'
    if 'property' in expected:
     props={x.get('name'):x.get('value') for x in testcase.findall('./properties/property')}
     assert props.get(expected['property']['name'])==expected['property']['value'],props
   else: assert testcase.get('name')=='OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims'
   row.update(accepted=True,owned_wall_seconds=owner['wall_seconds'],peak_sampled_rss_kib=owner['peak_sampled_rss_kib'],xml_case=testcase.get('name'))
  except Exception as error:
   row['review_error']=repr(error)
  rows.append(row)
  (work/'focused-result.json').write_text(json.dumps({'passed':len(rows)==6 and all(x['accepted'] for x in rows),'source_revision':'4498cc41168ad29fc4c111819a4cb7ccaa54a8bf','completed_case_count':len(rows),'required_case_count':6,'cases':rows,'qualification_claimed':False,'scope':'Finite exact-floor/shared adaptive and original Ellis control only; full science/frame/platform/release qualification remains incomplete.'},indent=2)+'\n')
  print(json.dumps(row),flush=True)
  if not row['accepted']:sys.exit(1)
print('SIX_FOCUSED_CONTROLS_COMPLETE')
