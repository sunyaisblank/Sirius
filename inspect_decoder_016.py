"""Inspect every 0.1.6 array and disclose differences from the earlier decoder."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import collections,hashlib,json,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(name):return json.loads((W/name).read_text())
n=read('isa016-result.json');old=read('isa-result.json');m=read('saved-isa-map.json')
inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),W/'isa016-result.json',W/'isa016-run.json',W/'isa-result.json',W/'saved-isa-map.json',W/'decoder016-types.h',W/'decoder016-release.json')}
lookup={s['base']+i['offset']:(s['stage'],i['instruction']) for s in m['stages'] for i in s['instructions']}
def array(v):
 data=bytearray()
 for entry in v['files']:
  p=ROOT/entry['path'];assert seal(p)=={k:entry[k] for k in ('bytes','sha256')};inputs[entry['path']]=seal(p);data.extend(p.read_bytes())
 assert len(data)==v['records']*v['record_bytes'];return data
def instructions(v):return list(struct.iter_unpack('<IiqQQ',array(v['arrays']['instructions'])))
waves=[]
for ordinal,wave in enumerate(n['wave_records']):
 inst=instructions(wave);prior=instructions(old['wave_records'][ordinal]);categories=collections.Counter();stages=collections.Counter();zero=0;unknown=0;negative=0;backwards=0;prev=None
 for cat,duration,time,addr,marker in inst:
  assert cat&255<13;categories[str(cat&255)]+=1;negative+=duration<0;backwards+=prev is not None and time<prev;prev=time
  if marker==0 and addr in lookup:stages[lookup[addr][0]]+=1
  elif marker==0 and addr==0:zero+=1
  else:unknown+=1
 timeline=list(struct.iter_unpack('<ii',array(wave['arrays']['timeline'])));assert all(0<=state<6 and duration>=0 for state,duration in timeline)
 total=sum(duration for state,duration in timeline)
 boundaries=[v for v in n['occupancy_records'] if (v['cu'],v['simd'],v['wave_id'])==(wave['cu'],wave['simd'],wave['wave_id']) and v['time'] in (wave['begin'],wave['end'])]
 assert any(v['start'] and v['time']==wave['begin'] for v in boundaries)
 assert any(not v['start'] and v['time']==wave['end'] for v in boundaries)
 assert total==wave['end']-wave['begin']
 waves.append({'ordinal':ordinal,**{k:wave[k] for k in ('cu','simd','wave_id','begin','end')},'instructions':len(inst),'categories':dict(categories),'bound_source_records':dict(stages),'zero_pc_records':zero,'other_unbound_records':unknown,'negative_duration_records':negative,'backward_time_records':backwards,'timeline_records':len(timeline),'timeline_sum':total,'declared_lifetime':wave['end']-wave['begin'],'occupancy_boundaries':boundaries,'prior_record_count_equal':len(inst)==len(prior),'prior_ordered_pc_values_equal':[(v[3],v[4]) for v in inst]==[(v[3],v[4]) for v in prior],'prior_category_differences':sum((a[0]&255)!=(b[0]&255) for a,b in zip(inst,prior))})
other=[]
for part in n['other_simd_records']:
 for size,time,cycles,wgp,cat,pad in struct.iter_unpack('<QqHBB4s',array(part)):
  assert size==24 and cat<13;other.append({'time':time,'execution_cycles':cycles,'encoded_wgp':wgp,'category':cat})
assert len(other)==2663
other_summary={'records':len(other),'first_issue_clock':min(v['time'] for v in other),'last_issue_clock':max(v['time'] for v in other),'encoded_wgp_counts':dict(collections.Counter(str(v['encoded_wgp']) for v in other)),'category_counts':dict(collections.Counter(str(v['category']) for v in other)),'sum_execution_cycles_inclusive':sum(v['execution_cycles'] for v in other),'max_execution_cycles':max(v['execution_cycles'] for v in other),'backward_issue_records':sum(a['time']>b['time'] for a,b in zip(other,other[1:])),'scope':'Public struct has no PC, wave or SIMD identifier. Cycles exclude stall; operations may overlap. This cannot identify exclusive shader/source/event cost.'}
assert set(v['encoded_wgp'] for v in other)=={0}
assert len(n['realtime_records'])==len(old['realtime_records'])==11912
rt_deltas=[a['shader_clock']-b['shader_clock'] for a,b in zip(n['realtime_records'],old['realtime_records'])]
rt={'records':len(rt_deltas),'ordered_reference_clock_values_equal':all(a['realtime_clock']==b['realtime_clock'] for a,b in zip(n['realtime_records'],old['realtime_records'])),'shader_clock_equal_records':sum(d==0 for d in rt_deltas),'shader_clock_delta_min':min(rt_deltas),'shader_clock_delta_max':max(rt_deltas),'distinct_shader_clock_deltas':len(set(rt_deltas)),'rt_frequency_records':n['rt_frequency_hz_records'],'scope':'Published gfx11 idle/stall correction changes decoded shader clocks. Earlier decoder clock-alignment observation remains historical and is not reusable for this result. No reference-clock frequency or exact UI/SPM/event identity established.'}
assert rt['ordered_reference_clock_values_equal']
occ=[]
assert len(n['occupancy_records'])==len(old['occupancy_records'])==72
for new,prior in zip(n['occupancy_records'],old['occupancy_records']):
 unchanged=[k for k in prior if k not in ('cu','time')];assert all(new[k]==prior[k] for k in unchanged)
 occ.append({'prior_encoded_cu':prior['cu'],'current_encoded_cu':new['cu'],'current_shader_array':(new['cu']&128)>>7,'current_wgp':new['cu']&127,'simd':new['simd'],'start':new['start'],'prior_time':prior['time'],'current_time':new['time'],'addr':new['addr']})
result={'source_revision':n['source_revision'],'input_seals':inputs,'status':n['status'],'info_records':n['info_records'],'isa_status_counts':n['isa_status_counts'],'unknown_requests':n['first_unknown_isa_addresses'],'waves':waves,'other_simd':other_summary,'realtime_comparison':rt,'occupancy_comparison':occ,'record_counts':n['record_counts'],'scope':'Actual full-array observations on unchanged trace/source. Absence of INFO does not establish complete stitching/source coverage: auxiliary ISA remains unbound with 16 zero PCs. Detailed target is still WGP0/SIMD0. New other-SIMD records are anonymous selected-WGP issue observations; no dominant instruction cost, supported hot-lane capture control, optimization or full qualification.'}
p=W/'decoder016-inspection.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'waves':[{k:v[k] for k in ('ordinal','instructions','bound_source_records','zero_pc_records','timeline_sum','declared_lifetime')} for v in waves],'other_simd':other_summary,'realtime_comparison':rt}))
