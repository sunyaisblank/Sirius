"""Bind owned saved-profile observations to the preserved code/source join."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from decimal import Decimal
import hashlib,json
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
paths=[Path(__file__),W/'trace-code-join.json',W/'counter-view-plan.json',W/'clock-domain-review.json']
join=read(paths[1]);plan=read(paths[2]);trace=ROOT/plan['protected_trace']['path'];paths.append(trace)
assert seal(trace)=={k:plan['protected_trace'][k] for k in ('bytes','sha256')}
def elements(path):
 paths.append(path);r=read(path);return r['elements'] if isinstance(r,dict) else r
def text(rows,suffix):
 values=[e['name'] for e in rows if e['control']=='ControlType.Text' and e['automation_id'].endswith('.'+suffix)]
 assert values and len(set(values))==1,(suffix,values)
 return values[0]
def micro(value):
 assert value.endswith(' μs');return Decimal(value[:-3].replace(',',''))
events=[]
for number,stage,expected_api in ((16,'Transport','b4bad5049f1d989e'),(22,'Endpoint','a08a2f262eca2b69')):
 action=f'profiler-counter-event{number}';bridge=W/(action+'-bridge.json');paths.append(bridge);r=read(bridge)
 assert r['accepted'] and r['source_revision']==join['source_revision']
 o=r['Windows_owner'];assert o['selected_event_confirmed'] and o['trace_sha256']==plan['trace']['sha256'] and o['owned_process_absent'] and not o['cleanup_errors']
 event=elements(W/action/'event-automation.json');pipeline=elements(W/action/'pipeline-automation.json');region=elements(W/action/'region-automation.json')
 api=text(event,'api_pso_hash_value_')[2:].lower();internal=text(event,'internal_pipeline_hash_value_')[2:].lower()
 assert api==expected_api
 matches=[v for v in join['pso_correlations'] if v['api_pso']==api and v['internal_upper']+v['internal_lower']==internal]
 assert len(matches)==1
 start,end,duration=(text(event,k) for k in ('start_time_','end_time_','duration_'))
 assert micro(end)-micro(start)==micro(duration)==micro(text(event,'work_duration_'))
 region_start,region_end,region_duration=(text(region,k) for k in ('start_time_','end_time_','duration_'))
 assert micro(region_end)-micro(region_start)==micro(region_duration) and micro(start)<micro(region_start)<micro(region_end)<micro(end)
 hashes={}
 for observation in o['observations']:
  paths.append(W/action/(observation['label']+'-automation.json'))
  images=observation['window_observation'];assert len(images)==1 and images[0]['print_window_succeeded']
  image=images[0];p=W/action/Path(image['image_path'].replace('\\','/')).name;paths.append(p)
  assert seal(p)=={'bytes':image['image_bytes'],'sha256':image['image_sha256']}
  hashes[observation['label']]=image['image_sha256']
 assert all(hashes[k]==hashes['counter-graphs'] for k in ('cache-tooltip','lds-tooltip','memory-tooltip'))
 static={k:text(pipeline,k) for k in ('total_wavefronts_','total_threads_','wavefront_mode_','vector_registers_','scalar_registers_','local_data_share_per_tg_','limiting_occupancy_','theoretical_occupancy_','resource_savings_for_extra_wave_')}
 events.append({'event':text(event,'event_name_'),'stage':stage,'api_pso':api,'internal_pipeline':internal,'code_object_ordinal':join[stage]['code_object']['ordinal'],'start':start,'end':end,'duration':duration,'explicit_region':{'start':region_start,'end':region_end,'duration':region_duration},'pipeline_static_observation':static,'counter_tooltips_observed':False,'scope':'RGP-rendered instrumented trace observation and exact PSO/internal-code join. Region/event selection does not establish exclusive counter ownership; theoretical occupancy is not executed bottleneck proof.'})
failed=W/'profiler-counter-hover16-bridge.json';paths.append(failed);r=read(failed);assert not r['accepted'] and r['returncode']==1 and not r['errors'] and r['linux_child_birth_absent']
owner=W/'profiler-counter-hover16/owner.json';paths.append(owner);o=read(owner)
assert o['status']=='failed' and 'Cursor position unavailable' in o['exception'] and 'cursor' not in o and not o['cleanup_errors'] and o['owned_process_absent']
recovery=W/'profiler-counter-hover16/cursor-recovery.json';paths.append(recovery);c=read(recovery);assert c['pass'] and c['result']=='no-cursor-token'
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
result={'source_revision':join['source_revision'],'input_seals':inputs,'trace':seal(trace),'events':events,'real_hover_attempt':{'accepted':False,'exception':o['exception'],'cursor_move_performed':False,'recovery':c,'next_attempt_rejected_without_changed_premise':True},'clock_domain':'SPM REFCLK timestamp counts, configured shader-cycle sampling interval, SQTT-rendered event times and Vulkan ordered-marker times have no established direct interval mapping. No nearest-sample or derived-value-delta attribution.','formula_scope':'Finite/cache-arithmetic checks are consistency observations; installed Phoenix/GPA formula and filter attestation remains open. PCIe bytes are a writer placeholder.','qualification':'No production source change, scientific dispatch, new trace acquisition, ordinary performance benefit, completed image or full scientific/release qualification in this batch.'}
for p,value in inputs.items():assert seal(ROOT/p)==value
p=W/'counter-findings.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'events_joined':len(events),'numeric_tooltips':False,'real_cursor_move':False,'production_change':False}))
