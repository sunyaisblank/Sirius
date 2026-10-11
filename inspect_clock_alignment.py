"""Record actual complete clock-phase observations without certifying an event epoch."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import bisect,collections,hashlib,json,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan=read(W/'plan.json');trace=ROOT/plan['trace']['path'];paths=[Path(__file__),W/'plan.json',trace,W/'spm-index-result.json',W/'raw-spm-result.json',W/'decoder-result.json',W/'decoder-run.json',W/'trace_decoder_types.h']
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths};assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
idx=read(W/'spm-index-result.json');d=read(W/'decoder-result.json');run=read(W/'decoder-run.json');assert d['status']==0 and not d['callback_errors']
assert run['returncode']==0 and not run['errors'] and run['final_owned_group']==[] and not run['limit_reason']
c=next(v for v in idx['chunks'] if v['name']=='SpmSession');b=trace.read_bytes();raw=b[c['data_offset']:c['data_offset']+c['data_bytes']]
assert hashlib.sha256(raw).hexdigest()==c['data_sha256'];t=struct.unpack('<11911Q',raw)
pairs=d['realtime_records'];assert len(pairs)==11912 and all(v['reserved']==0 for v in pairs)
shader=[v['shader_clock'] for v in pairs];rt=[v['realtime_clock'] for v in pairs]
assert all(a<z for a,z in zip(shader,shader[1:])) and all(a<z for a,z in zip(rt,rt[1:])) and all(a<z for a,z in zip(t,t[1:]))
offset=18691697672192;assert offset==272*(1<<36)==17*(1<<40)
residuals=[];anchors=[];ties=[]
for sample,x in enumerate(t):
 u=x-offset;i=bisect.bisect_left(rt,u);candidates=list(range(max(0,i-1),min(len(rt),i+1)));assert candidates
 distances=[abs(rt[z]-u) for z in candidates];best=min(distances);nearest=[z for z in candidates if abs(rt[z]-u)==best]
 if len(nearest)>1:ties.append(sample)
 z=nearest[0];residuals.append(u-rt[z]);anchors.append(z)
occupancy=[]
for addr in sorted({v['addr'] for v in d['occupancy_records']}):
 rows=[v for v in d['occupancy_records'] if v['addr']==addr]
 occupancy.append({'raw_address_hex':hex(addr),'record_count':len(rows),'start_records':sum(v['start'] for v in rows),'end_records':sum(not v['start'] for v in rows),'first_shader_time':min(v['time'] for v in rows),'last_shader_time':max(v['time'] for v in rows),'marker_ids':sorted({v['marker_id'] for v in rows}),'observed_cu_ids':sorted({v['cu'] for v in rows})})
result={'source_revision':plan['source_revision'],'input_seals':inputs,'trace':seal(trace),'status':d['status'],'info_records':d['info_records'],'info_strings':d['info_strings'],'clock_pairs':len(pairs),'shader_clock_strictly_increasing':True,'realtime_clock_strictly_increasing':True,'spm_clock_strictly_increasing':True,'rt_frequency_hz_records':d['rt_frequency_hz_records'],'first_pair':pairs[0],'last_pair':pairs[-1],'paired_shader_delta':shader[-1]-shader[0],'paired_reference_delta':rt[-1]-rt[0],'whole_span_shader_units_per_reference_unit':(shader[-1]-shader[0])/(rt[-1]-rt[0]),'candidate_reference_epoch_offset':offset,'offset_ambiguity':'Same observed offset is 272 * 2^36 and 17 * 2^40; these bytes alone do not identify the truncated timestamp width.','nearest_anchor_residual_counts':dict(sorted(collections.Counter(residuals).items())),'maximum_absolute_reference_residual':max(map(abs,residuals)),'distinct_nearest_anchors':len(set(anchors)),'nearest_anchor_tie_samples':ties,'nearest_anchor_indices_strictly_increasing':all(a<z for a,z in zip(anchors,anchors[1:])),'nonzero_or_minus4_residuals':[{'sample':i,'reference_residual':v,'nearest_anchor':anchors[i]} for i,v in enumerate(residuals) if v not in (0,-4)],'occupancy_address_groups':occupancy,'scope':'Full observed finite clock-phase proximity only. No exact sample identity, independently certified REFCLK epoch/frequency binding, RGP displayed-event origin, complete ISA stitching, exclusive counter ownership, ordinary benefit or scientific/release qualification.'}
for p,v in inputs.items():assert seal(ROOT/p)==v
p=W/'clock-alignment-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ('clock_pairs','paired_shader_delta','paired_reference_delta','candidate_reference_epoch_offset','nearest_anchor_residual_counts','maximum_absolute_reference_residual','distinct_nearest_anchors','nearest_anchor_indices_strictly_increasing')}))
