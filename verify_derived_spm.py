"""Read documented derived SPMv1 shapes, identities and component/group links."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, math, struct, subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*a):return subprocess.check_output(['git',*a],text=True).strip()
head=read(WORK/'source.json')['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
refs=read(WORK/'derived-references.json')
for v in refs['files']:assert seal(ROOT/v['path'])=={k:v[k] for k in ('bytes','sha256')}
index=read(WORK/'spm-index-result.json');raw=read(WORK/'raw-spm-result.json')
assert raw['trace']==index['trace'] and raw['session']['strictly_increasing']
trace=WORK/'native-spm-capture-console/spm-original-legacy-interval.rgp';assert seal(trace)==index['trace'];b=trace.read_bytes()
paths=[Path(__file__),trace,WORK/'spm-index-result.json',WORK/'raw-spm-result.json',WORK/'derived-references.json']+[ROOT/v['path'] for v in refs['files']]
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
n=raw['session']['timestamps'];sessions=[c for c in index['chunks'] if c['name']=='DerivedSpmSes'];assert len(sessions)==1
s=sessions[0];assert s['version']==1 and s['compression']==0 and s['header_bytes']==0 and s['data_bytes']==12
pci,expected_counters,expected_groups=struct.unpack_from('<3I',b,s['data_offset']);assert (pci,expected_counters,expected_groups)==(0,33,5)
counters={};groups={}
def text(data):
 assert b'\0' not in data
 return data.decode('utf-8')
for c in index['chunks']:
 if c['name'] not in ('DerivedSpmCtr','DerivedSpmGroup'):continue
 assert c['version']==1 and c['compression']==0
 data=b[c['data_offset']:c['data_offset']+c['data_bytes']]
 if c['name']=='DerivedSpmCtr':
  assert c['header_bytes']==24
  cpci,idx,components,unit,nl,dl=struct.unpack_from('<6I',b,c['header_offset'])
  assert cpci==pci and idx<expected_counters and idx not in counters and components<=expected_counters and nl>0
  assert len(data)==n*8+components*4+nl+dl
  values=struct.unpack_from('<'+str(n)+'d',data)
  offset=n*8;links=list(struct.unpack_from('<'+str(components)+'I',data,offset));offset+=components*4
  assert len(set(links))==len(links) and idx not in links and all(k<expected_counters for k in links)
  name=text(data[offset:offset+nl]);offset+=nl;description=text(data[offset:offset+dl])
  finite=[v for v in values if math.isfinite(v)]
  counters[idx]={'counter_index':idx,'chunk_ordinal':c['ordinal'],'display_name':name,'description':description,
                 'usage_type_encoding':unit,'component_indices':links,'samples':n,
                 'sample_data_sha256':hashlib.sha256(data[:n*8]).hexdigest(),
                 'finite_samples':len(finite),'nonfinite_sample_indices':[j for j,v in enumerate(values) if not math.isfinite(v)],
                 'nonzero_finite_samples':sum(v!=0 for v in finite),
                 'minimum_finite':min(finite) if finite else None,'maximum_finite':max(finite) if finite else None}
 else:
  assert c['header_bytes']==20
  gpci,idx,members,nl,dl=struct.unpack_from('<5I',b,c['header_offset'])
  assert gpci==pci and idx<expected_groups and idx not in groups and 0<members<=expected_counters and nl>0
  assert len(data)==members*4+nl+dl
  links=list(struct.unpack_from('<'+str(members)+'I',data));offset=members*4
  assert len(set(links))==len(links) and all(k<expected_counters for k in links)
  name=text(data[offset:offset+nl]);offset+=nl;description=text(data[offset:offset+dl])
  groups[idx]={'group_index':idx,'chunk_ordinal':c['ordinal'],'name':name,'description':description,'member_indices':links}
assert set(counters)==set(range(expected_counters)) and set(groups)==set(range(expected_groups))
for c in counters.values():c['component_display_names']=[counters[k]['display_name'] for k in c['component_indices']]
for g in groups.values():g['member_display_names']=[counters[k]['display_name'] for k in g['member_indices']]
assert {g['name'] for g in groups.values()}=={'Cache','Ray tracing','LDS','Memory (%)','Memory (bytes)'}
for p,v in inputs.items():assert seal(ROOT/p)==v
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
result={'source_revision':head,'input_seals':inputs,'trace':index['trace'],
        'session':{'pci_id':pci,'raw_session_pci_id':raw['session']['pci_id'],'pci_id_scope':'Public writer deliberately writes zero for derived chunks; do not equate this placeholder with raw/device PCI identity.','counter_count':expected_counters,'group_count':expected_groups,'samples_per_counter':n},
        'counters':list(counters.values()),'groups':list(groups.values()),
        'scope':'Recorded tool-derived binary64 values, length-counted UTF8 strings, counter_index identities and complete component/group links conform to the pinned public v1 writer schema. No formula reconstruction or unit/interval semantic verification, event/pipeline/source cost join, exclusive bottleneck, hardware qualification, ordinary performance benefit, rendered frame or full release/scientific qualification claim.'}
p=WORK/'derived-spm-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'counter_count':len(counters),'group_count':len(groups),'samples_per_counter':n,
                  'finite_counter_streams':sum(c['finite_samples']==n for c in counters.values()),
                  'groups':[{k:g[k] for k in ('name','member_display_names')} for g in groups.values()]}))
