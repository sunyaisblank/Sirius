"""Independent arithmetic on original raw XML; never execute/resample a workload."""
from pathlib import Path
from decimal import Decimal,localcontext
import json,xml.etree.ElementTree as ET
R=Path.cwd();W=Path(__file__).resolve().parent;RUN=R/'out/native-original-input-ownership/runtime'
def props(label):
 root=ET.parse(RUN/label/'gtest.xml').getroot()
 assert root.get('tests')=='1' and root.get('failures')==root.get('errors')=='0'
 assert not any(root.findall('.//'+tag) for tag in ['failure','error','skipped'])
 out={}
 for row in root.findall('.//property'):
  k=row.get('name');assert k not in out;out[k]=row.get('value')
 return out
raw={label:props(label) for label in ['a1','b1','b2','a2']}
accepted=json.loads((W/'benefit.json').read_text());assert accepted['audit_passed']
expected={(x['direction'],x.get('family','instance'),x['scope']):x['pass'] for x in accepted['checks']+accepted['accounting_checks']}
checks=[]
def number(s):
 v=Decimal(s);assert v.is_finite() and v>=0;return v
def median(v):assert len(v)==3;return sorted(v)[1]
def check(direction,family,scope,a,b,minimum):
 with localcontext() as context:
  context.prec=60
  gain=100*(a-b)/a;passed=gain>=minimum
  assert expected[(direction,family,scope)]==passed
  checks.append({'direction':direction,'family':family,'scope':scope,'baseline_ms':str(a),'candidate_ms':str(b),'improvement_percent':str(gain),'minimum_percent':str(minimum),'passed':passed})
for direction,al,bl in [('forward','a1','b1'),('reverse','a2','b2')]:
 for family,prefix,stage_count in [('legacy','legacy_repeat',3),('modern','coupled_repeat',5)]:
  def whole(p):return median([number(p[f'{prefix}{r}_complete_call_ms']) for r in range(3)])
  check(direction,family,'complete_call_ms',whole(raw[al]),whole(raw[bl]),Decimal(5))
  def stages(p):
   totals={}
   for i in range(51):
    fields=dict(item.split('=',1) for item in p[f'device_observation_{i}'].split(';'))
    if not fields['phase'].startswith(prefix):continue
    r=int(fields['phase'][len(prefix):]);assert 0<=r<3
    totals.setdefault(fields['stage'],[Decimal(0)]*3)[r]+=number(fields['total_ms'])
   assert len(totals)==stage_count
   return totals
  a,b=stages(raw[al]),stages(raw[bl]);assert a.keys()==b.keys()
  for stage in sorted(a):check(direction,family,stage+' inclusive Dispatch',median(a[stage]),median(b[stage]),Decimal(-5))
 def api(p):
  keys=['observed_create_ms','disabled_attempt_ms','disabled_sample_ms']
  keys += [f'{prefix}{r}_complete_call_ms' for prefix in ['legacy_repeat','coupled_repeat'] for r in range(3)]
  assert p['owner_provisioning_ms']=='0'
  return sum((number(p[k]) for k in keys),Decimal(0))
 check(direction,'instance','observed-instance Create and all11 original numerical API windows',api(raw[al]),api(raw[bl]),Decimal(0))
assert len(checks)==22
passed=all(x['passed'] for x in checks);assert passed==accepted['benefit_passed']
out=W/'decimal-agreement.json';assert not out.exists()
out.write_text(json.dumps({'agrees_on_all22':True,'all22_pass':passed,'checks':checks,'scope':'Independent Decimal/raw-XML numerical witness only. No additional execution, provenance or scientific acceptance.'},indent=2)+'\n')
print(json.dumps({'agrees_on_all22':True,'checks_passed':sum(x['passed'] for x in checks),'checks_total':22}))
