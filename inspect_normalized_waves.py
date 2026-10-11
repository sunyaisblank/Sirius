"""Inspect every retained instruction/timeline record; keep partial coverage explicit."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import collections,hashlib,json,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
r=json.loads((W/'isa2-result.json').read_text());plan=json.loads((W/'plan.json').read_text());prior=ROOT/plan['prior_decoder_result']['path'];assert seal(prior)=={k:plan['prior_decoder_result'][k] for k in ('bytes','sha256')};old=json.loads(prior.read_text())
inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),W/'isa2-result.json',W/'isa2-run.json',W/'saved-isa2-map.json',W/'trace_decoder_types.h',prior)}
m=json.loads((W/'saved-isa2-map.json').read_text());lookup={v['base']+i['offset']:{'stage':v['stage'],**i} for v in m['stages'] for i in v['instructions']}
def arrays(v):
 b=bytearray()
 for entry in v['files']:
  p=ROOT/entry['path'];assert seal(p)=={k:entry[k] for k in ('bytes','sha256')};inputs[entry['path']]=seal(p);b.extend(p.read_bytes())
 assert len(b)==v['records']*v['record_bytes'];return b
waves=[]
for v in r['wave_records']:
 b=arrays(v['arrays']['instructions']);categories=collections.Counter();resolved=0;zero=0;unknown=[];bound=[];first_unbound=None;previous=None;backwards=0;negative=0
 for ordinal,offset in enumerate(range(0,len(b),32)):
  cat,duration,clock,addr,marker=struct.unpack_from('<IiqQQ',b,offset);category=cat&255;stall=cat>>8;assert category<13;categories[str(category)]+=1
  item={'ordinal':ordinal,'category':category,'stall':stall,'duration':duration,'time':clock,'addr':addr,'marker_id':marker}
  if previous is not None and clock<previous:backwards+=1
  previous=clock;negative+=duration<0
  if marker==0 and addr in lookup:
   resolved+=1
   if len(bound)<32:bound.append({**item,'saved_instruction':lookup[addr]})
  elif addr==0 and marker==0:
   zero+=1
   if first_unbound is None:first_unbound=item
  elif len(unknown)<32:unknown.append(item)
 t=arrays(v['arrays']['timeline']);states=collections.Counter();duration=0
 for offset in range(0,len(t),8):
  state,cycles=struct.unpack_from('<ii',t,offset);assert 0<=state<6;states[str(state)]+=1;duration+=cycles
 wave={k:v[k] for k in ('cu','simd','wave_id','contexts','begin','end')};wave.update(instructions=v['arrays']['instructions']['records'],categories=dict(categories),resolved_pc_records=resolved,zero_pc_records=zero,other_unresolved_examples=unknown,first_zero_pc_record=first_unbound,bound_examples=bound,backward_instruction_times=backwards,negative_duration_records=negative,timeline_states=dict(states),timeline_duration_sum=duration,declared_wave_duration=v['end']-v['begin'])
 waves.append(wave)
unchanged={k:r[k]==old[k] for k in ('realtime_records','occupancy_records','rt_frequency_hz_records')};assert all(unchanged.values())
result={'source_revision':r['source_revision'],'input_seals':inputs,'status':r['status'],'isa_status_counts':r['isa_status_counts'],'info_strings':r['info_strings'],'prior_records_whole_value_equal':unchanged,'prior_record_counts':old['record_counts'],'current_record_counts':r['record_counts'],'wave_array_bytes':r['wave_array_bytes'],'waves':waves,'scope':'Actual partial selected-wave observation; zero-PC/category timing cannot be charged to source lines or UI events. Incomplete stitching, missing frequency and display origin remain. No exclusive cost, bottleneck, optimized code or qualification.'}
p=W/'normalized-wave-inspection.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'unchanged':unchanged,'waves':[(v['instructions'],v['resolved_pc_records'],v['zero_pc_records'],v['first_zero_pc_record']) for v in waves]}))
