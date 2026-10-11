"""Validate documented raw SPM structures and coverage; do not attribute event cost."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
from collections import Counter, defaultdict
import hashlib, json, struct, subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*a):return subprocess.check_output(['git',*a],text=True).strip()
head=read(WORK/'source.json')['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
trace=WORK/'native-spm-capture-console/spm-original-legacy-interval.rgp'
index=read(WORK/'spm-index-result.json');assert index['source_revision']==head and seal(trace)==index['trace']
refs=read(WORK/'references.json')
for r in refs['files']:assert seal(ROOT/r['path'])=={k:r[k] for k in ('bytes','sha256')}
paths=[Path(__file__),trace,WORK/'spm-index-result.json',WORK/'references.json',WORK/'pal-gpu-perf-header.h',WORK/'pal-gpu-perf-producer.cpp']
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
b=trace.read_bytes();sessions=[c for c in index['chunks'] if c['name']=='SpmSession'];assert len(sessions)==1
s=sessions[0];assert s['version']==2 and s['compression']==0 and s['header_bytes']==20
pci,flags,freq,n,expected=struct.unpack_from('<5I',b,s['header_offset'])
assert flags==0 and freq==index['spm_configuration']['sampleFrequency']==4096 and n>1 and s['data_bytes']==n*8
timestamps=struct.unpack_from('<'+str(n)+'Q',b,s['data_offset'])
assert all(a<b for a,b in zip(timestamps,timestamps[1:]))
deltas=[b-a for a,b in zip(timestamps,timestamps[1:])]
requested={tuple(v) for v in index['spm_configuration']['perfCounters']}
assert index['spm_enabled'] and len(requested)==len(index['spm_configuration']['perfCounters'])==27
assert all(instance==4095 for _,instance,_ in requested)
rows=[];ids=set();actual=defaultdict(set)
for c in index['chunks']:
 if c['name']!='SpmCounterData':continue
 assert c['version']==2 and c['compression']==0 and c['header_bytes']==20
 cpci,block,instance,event,width=struct.unpack_from('<5I',b,c['header_offset'])
 assert cpci==pci and width in (2,4) and c['data_bytes']==n*width
 identity=(block,instance,event);assert identity not in ids;ids.add(identity)
 assert (block,4095,event) in requested
 actual[(block,event)].add(instance)
 fmt='H' if width==2 else 'I';values=struct.unpack_from('<'+str(n)+fmt,b,c['data_offset'])
 nonzero=sum(v!=0 for v in values)
 rows.append({'chunk_ordinal':c['ordinal'],'block':block,'instance':instance,'event':event,
              'item_bytes':width,'samples':n,'data_sha256':c['data_sha256'],
              'nonzero_samples':nonzero,'minimum':min(values),'maximum':max(values),
              'sum_raw_values':sum(values),'first_nonzero_sample':next((i for i,v in enumerate(values) if v),None),
              'last_nonzero_sample':next((i for i in range(n-1,-1,-1) if values[i]),None)})
assert len(rows)==expected==183
assert set(actual)=={(block,event) for block,_,event in requested}
for key,instances in actual.items():
 assert instances==set(range(len(instances))), (key,instances)
assert any(r['nonzero_samples'] for r in rows)
for p,v in inputs.items():assert seal(ROOT/p)==v
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
result={'source_revision':head,'input_seals':inputs,'trace':index['trace'],
        'session':{'pci_id':pci,'flags':flags,'sample_frequency_gpu_cycles':freq,'timestamps':n,'raw_counter_streams':expected,
                   'first_timestamp':timestamps[0],'last_timestamp':timestamps[-1],
                   'strictly_increasing':True,'timestamp_delta_counts':dict(sorted(Counter(deltas).items())),
                   'timestamp_data_sha256':s['data_sha256']},
        'requested_rows':len(requested),'expanded_counter_rows':rows,
        'requested_instance_coverage':[{'block':block,'event':event,'instances':sorted(instances)} for (block,event),instances in sorted(actual.items())],
        'raw_counter_streams_with_nonzero_values':sum(r['nonzero_samples']>0 for r in rows),
        'scope':'Documented raw SPMv2 size/type, all requested block/event-to-instance rows, nonempty strictly increasing timestamp coverage and raw-value summaries. Values have no assigned counter names/units or event/source-window meaning here. No exclusive eligibility cause, hardware qualification, bottleneck/optimization, ordinary timing benefit, rendered frame or release/scientific qualification claim.'}
p=WORK/'raw-spm-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'timestamps':n,'raw_counter_streams':len(rows),'requested_counter_rows':len(requested),
                  'nonzero_counter_streams':result['raw_counter_streams_with_nonzero_values'],
                  'timestamp_delta_counts':result['session']['timestamp_delta_counts']}))
