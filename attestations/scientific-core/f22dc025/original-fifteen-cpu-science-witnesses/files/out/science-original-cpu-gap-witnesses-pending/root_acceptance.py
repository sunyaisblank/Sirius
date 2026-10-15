from pathlib import Path
import hashlib,json,sys,subprocess,xml.etree.ElementTree as ET
r=Path.cwd();w=r/'out/science-original-cpu-gap-witnesses-pending';label=sys.argv[1]
def doc(p):return json.loads(p.read_text())
def pin(p):
 raw=p.read_bytes();return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
def check(p,v):
 assert pin(p)=={k:v[k] for k in ['bytes','sha256']},str(p)
 if 'resolved_path' in v:assert str(p.resolve())==v['resolved_path']
def absent(pid,ticks):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]!=str(ticks)
 except FileNotFoundError:return True
a=doc(w/'actions.json');order=list(a);assert label in order
for old in order[:order.index(label)]:
 prior=doc(w/('root-accepted-'+old+'.json'));assert prior['passed']
 for rel,v in prior['files'].items():check(r/rel,v)
pf=doc(w/'root-preflight.json')
for rel,v in pf['files'].items():check(Path(rel) if rel.startswith('/') else r/rel,v)
o=doc(w/label/'owner.json');assert o['passed'] is True and o['status']=='terminal' and o['child_exit']==0 and o['stop_reason'] is None and not o['cleanup_errors'] and o['remaining_owned_processes']==[]
assert o['source']=={'revision':a[label]['revision'],'status':''} and o['source_unchanged'] and o['whole_inputs_unchanged']
before=doc(w/label/'inputs-before.json');after=doc(w/label/'inputs-after.json');assert before==after
for rel,v in before.items():check(Path(rel) if rel.startswith('/') else r/rel,v)
assert o['observed_births'] and all(x['absent'] and absent(x['pid'],x['saved_start_ticks']) for x in o['observed_births'])
assert absent(o['controller']['pid'],o['controller']['start_ticks'])
if a[label].get('software_provider', False):
    assert o['provider_mappings'] and all(x['pid']==o['child_pid'] and str(x['start_ticks'])==str(o['child_start_ticks']) and x['owned_session']==o['child_pid'] and x['maps'] for x in o['provider_mappings'])
else:
    assert o['provider_mappings']==[]
assert o['peak_sampled_rss_kib']<=4194304 and o['minimum_sampled_available_kib']>=2097152 and o['wall_seconds']<300
x=ET.parse(w/label/'gtest.xml').getroot();cases=x.findall('.//testcase');e=a[label]['root_required_xml'];assert len(cases)==int(x.get('tests','-1'))==1
assert all(int(n.get(k,'0'))==0 for n in [x,*x.findall('.//testsuite')] for k in ['failures','errors','skipped','disabled'])
assert not any(x.findall('.//'+k) for k in ['failure','error','skipped'])
c=cases[0];assert c.get('classname')==e['suite'] and c.get('name')==e['name'] and c.get('status')=='run' and c.get('result')=='completed'
properties={n.get('name'):n.get('value') for n in c.findall('./properties/property')}
if 'property' in e:
    assert properties.get(e['property']['name'])==e['property']['value']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==a[label]['revision'] and not subprocess.check_output(['git','status','--porcelain']).strip()
paths=sorted(p for p in (w/label).rglob('*') if p.is_file());files={str(p.relative_to(r)):pin(p) for p in paths}
result={'schema':'sirius-kerr-shared-loop-root-accepted-control-v1','passed':True,'case':label,'source_revision':a[label]['revision'],'exact_one_completed_XML_zero_FESD':True,'required_positive_property':e.get('property'),'test_seconds':c.get('time'),'owner_seconds':o['wall_seconds'],'peak_sampled_aggregate_RSS_KiB':o['peak_sampled_rss_kib'],'fresh_whole_input_rows':len(before),'owned_child_and_controller_births_absent':True,'source_and_inputs_unchanged':True,'files':files,'scope':'One exact focused parity control with original assertions and guards; no fullimage/allP1/Mandatory/latency/release qualification.'}
p=w/('root-accepted-'+label+'.json');assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['case','passed','test_seconds','owner_seconds','peak_sampled_aggregate_RSS_KiB','fresh_whole_input_rows']}))

