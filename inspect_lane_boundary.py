"""Read actual saved target/config/load bytes and join only recorded occupancy."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import collections,hashlib,json,re,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan=json.loads((W/'plan.json').read_text());trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
inputs={str(p.relative_to(ROOT)):seal(p) for p in [Path(__file__),W/'plan.json',trace,W/'spm-index-result.json',W/'isa-result.json',W/'saved-isa-summary.json',W/'trace_decoder_types.h',W/'lane-source-preparation.json',*sorted(W.glob('lane-pal-*')),*sorted(W.glob('lane-rdp-*'))]}
index=json.loads((W/'spm-index-result.json').read_text());data=trace.read_bytes()
def chunk(name,part='data'):
 c=next(v for v in index['chunks'] if v['name']==name);assert c['compression']==0
 b=data[c[part+'_offset']:c[part+'_offset']+c[part+'_bytes']];assert len(b)==c[part+'_bytes'] and hashlib.sha256(b).hexdigest()==c[part+'_sha256'];return c,b
a,asic=chunk('AsicInfo');assert a['version']==3 and len(asic)==608
fields=[]
for name,offset,fmt,value in [('shaderEngines',64,'<i',1),('computeUnitPerShaderEngine',68,'<i',12),('simdPerComputeUnit',72,'<i',2)]:
 assert struct.unpack_from(fmt,asic,offset)[0]==value
 fields.append({'name':name,'relative_offset':offset,'file_offset':a['data_offset']+offset,'raw_hex':asic[offset:offset+4].hex(),'value':value})
assert struct.unpack_from('<HH',asic,444)==(63,63) and asic[448:572]==bytes(124)
fields.append({'name':'cuMask[0][0..1]','relative_offset':444,'file_offset':a['data_offset']+444,'raw_hex':asic[444:448].hex(),'value':[63,63]})
s,header=chunk('SqttData','header');assert s['version']==5 and len(header)==40
pci,se,version,spec,api,wgp,buffer,flags=struct.unpack_from('<6IQI',header)
assert se==0 and wgp==0 and flags&3==1
t,cfgbytes=chunk('TraceConfig');config=json.loads(cfgbytes.rstrip(b'\0'));cfg=next(v for v in config['sources'] if v['name']=='gpuperfexp')['config']['sqtt']
assert cfg=={'enabled':True,'memoryLimitInMb':75,'enableInstructionTokens':True,'enableExecPopTokens':False,'seMask':1}
co,loads=chunk('COLoadEvent');count=struct.unpack('<I',chunk('COLoadEvent','header')[1])[0];assert count==3 and len(loads)==count*40
loaded=[]
for off in range(0,len(loads),40):
 pci_id,event,base,lo,hi,stamp=struct.unpack_from('<II4Q',loads,off);assert event==0
 loaded.append({'pci_id':pci_id,'event_type':event,'base_address':base,'internal_lower':f'{lo:016x}','internal_upper':f'{hi:016x}','cpu_timestamp':stamp})
assert 0x10001b800<min(v['base_address'] for v in loaded)
summary=json.loads((W/'saved-isa-summary.json').read_text());bases={v['base']:v['stage'] for v in summary['stages']}
r=json.loads((W/'isa-result.json').read_text());starts=[{**v,'stage':bases.get(v['addr'],'unbound_auxiliary' if v['addr']==0x10001b800 else 'unbound')} for v in r['occupancy_records'] if v['start']]
target=[v for v in starts if (v['cu'],v['simd'])==(0,0)]
assert [(v['stage'],v['time']) for v in target]==[('Transport',32204),('unbound_auxiliary',12911339)]
assert not any(v['stage']=='Endpoint' for v in target)
transport=[v for v in starts if v['stage']=='Transport' and 8415365<=v['time']<=8415369]
assert [(v['cu'],v['simd']) for v in transport]==[(8,3),(9,1),(10,3)]
producer=(W/'lane-pal-gpuPerfExperimentTraceSource.cpp').read_text();section=producer[producer.index('void GpuPerfExperimentTraceSource::OnSqttConfigUpdated('):producer.index('void GpuPerfExperimentTraceSource::OnSpmConfigUpdated(')]
keys=re.findall(r'GetValueByKey\("([^"]+)"',section);assert keys==['enabled','memoryLimitInMb','enableInstructionTokens','enableExecPopTokens','seMask']
gpa=(W/'lane-pal-gpaSession.cpp').read_text();region=gpa[gpa.index('ThreadTraceInfo sqttInfo = { };'):gpa.index('ThreadTraceTokenConfig tokenConfig =',gpa.index('ThreadTraceInfo sqttInfo = { };'))]
assert all(x not in region for x in ('threadTraceTargetSh','threadTraceTargetCu','threadTraceSimdMask'))
pal=(W/'lane-pal-gfx9PerfExperiment.cpp').read_text();assert '? traceInfo.optionValues.threadTraceTargetSh : 0' in pal and '? traceInfo.optionValues.threadTraceSimdMask : 0' in pal
result={'source_revision':plan['source_revision'],'input_seals':inputs,'asic_info':{'chunk':a,'fields':fields,'remaining_cu_masks_zero':True},'sqtt_header':{'chunk':s,'pci_id':pci,'shader_engine':se,'version':version,'instrumentation_spec':spec,'instrumentation_api':api,'reported_wgp_index':wgp,'trace_buffer_bytes':buffer,'instruction_timing_enabled':True,'exec_pop_tokens_enabled':False},'recorded_sqtt_config':cfg,'recorded_loads':loaded,'auxiliary_below_all_recorded_load_bases':True,'all_occupancy_starts':starts,'selected_target_starts':target,'second_transport_cluster':transport,'accepted_sqtt_parser_keys':keys,'supported_control_disposition':'Pinned RDP API/CLI/serializer and PAL gpuperfexp route provide no supported target-SA/CU/SIMD selector. Lower PAL AddThreadTrace supports them, but this GpaSession path leaves them unset. Private DefaultSqttDetailedCuIndex is not an established per-owned-RDP SIMD/SA control. No guessed keys or numeric mapping accepted.','scope':'Actual byte and occupancy observation plus bounded primary-source inspection. Raw 0.1.3 CU encoding retained. Cluster-to-UI-event identity remains an ordered-trial inference; no exclusive cost or optimization. Auxiliary ISA cannot create detailed timing for omitted targets.'}
p=W/'lane-boundary-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'starts':len(starts),'target_starts':[(v['stage'],v['time']) for v in target],'no_endpoint_detailed_target':True,'sqtt_keys':keys,'auxiliary_bound':False}))
