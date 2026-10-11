"""Inspect the exact acquired RDF index/configuration without decoding counter values."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, struct, subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*a):return subprocess.check_output(['git',*a],text=True).strip()
head=read(WORK/'source.json')['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
bridge=read(WORK/'native-spm-capture-console-bridge.json');assert bridge['accepted'] and not bridge['errors']
for p,v in bridge['input_seals'].items():assert seal(Path(p) if Path(p).is_absolute() else ROOT/p)==v,p
output=WORK/'native-spm-capture-console';o=read(output/'owner.json')
from verify_spm_capture_control import verify
control=verify(output)
trace=output/'spm-original-legacy-interval.rgp';identity=seal(trace)
assert identity=={k:o['trace'][k] for k in ('bytes','sha256')}
inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),trace,output/'owner.json',WORK/'native-spm-capture-console-bridge.json',WORK/'verify_spm_capture_control.py',WORK/'capture-spm-plan.json',WORK/'references.json')}
b=trace.read_bytes();magic,version,reserved,index,size=struct.unpack_from('<8sIIqq',b)
assert magic==b'AMD_RDF ' and version==3 and reserved==0
assert index>=32 and size>0 and size%64==0 and index+size==len(b) and size//64<=4096
chunks=[];spans=[];configs=[];infos=[]
for ordinal,start in enumerate(range(index,index+size,64)):
 name,compression,padding,v,ho,hs,do,ds,us=struct.unpack_from('<16sB3sIqqqqq',b,start)
 assert padding==b'\0'*3 and compression in (0,1) and v>0
 trimmed=name.rstrip(b'\0');assert trimmed and b'\0' not in trimmed;name=trimmed.decode('utf-8')
 for off,length in ((ho,hs),(do,ds)):
  assert 0<=off<=index and 0<=length<=index-off
  if length:assert off>=32;spans.append((off,off+length))
 if compression==0:assert us==0
 else:assert 0<us<=268435456
 h=b[ho:ho+hs];data=b[do:do+ds]
 c={'ordinal':ordinal,'name':name,'version':v,'compression':compression,'header_offset':ho,'header_bytes':hs,'header_sha256':hashlib.sha256(h).hexdigest(),'data_offset':do,'data_bytes':ds,'data_sha256':hashlib.sha256(data).hexdigest(),'uncompressed_size':us}
 if name=='TraceConfig':
  assert v==1 and hs==0 and compression==0 and data[-1:]==b'\0' and b'\0' not in data[:-1]
  config=json.loads(data[:-1].decode('utf-8'));configs.append(config);c['config']=config;c['terminal_nul_bytes']=1
 if name=='SystemInfo':
  assert v==1 and hs==0 and compression==0;info=json.loads(data.decode('utf-8'));infos.append(info);c['system_info']=info
 if name in ('SpmSession','SpmCounterData'):
  assert v==2 and hs==20;c['header_u32_le']=list(struct.unpack('<5I',h))
 chunks.append(c)
assert len(configs)==len(infos)==1
spans.sort();cursor=32;gaps=[]
for start,end in spans:
 assert start>=cursor and end>start
 if start>cursor:gaps.append({'offset':cursor,'bytes':start-cursor,'sha256':hashlib.sha256(b[cursor:start]).hexdigest()})
 cursor=end
if cursor<index:gaps.append({'offset':cursor,'bytes':index-cursor,'sha256':hashlib.sha256(b[cursor:index]).hexdigest()})
perf=[v for v in configs[0]['sources'] if v['name']=='gpuperfexp'];assert len(perf)==1
spm=perf[0]['config']['spm'];assert isinstance(spm['enabled'],bool)
for p,v in inputs.items():assert seal(ROOT/p)==v
result={'source_revision':head,'input_seals':inputs,'original_control':control,'trace':identity,'index_offset':index,'index_bytes':size,'chunks':chunks,'unindexed_regions':gaps,'spm_configuration':spm,'spm_enabled':spm['enabled'],'raw_spm_chunk_counts':{n:sum(c['name']==n for c in chunks) for n in ('SpmSession','SpmCounterData','SpmDb','DerivedSpmDb')},'scope':'Exact whole trace/index, recorded counter configuration and raw SPM headers. No counter value/schema/semantic or event/source cost attribution yet; no exclusive eligibility cause, hardware-support, ordinary performance, frame or qualification claim.'}
p=WORK/'spm-index-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'trace':identity,'indexed_chunks':len(chunks),'spm_enabled':result['spm_enabled'],'requested_counter_rows':len(spm.get('perfCounters',[])),'raw_spm_chunk_counts':result['raw_spm_chunk_counts'],'unindexed_regions':len(gaps),'chunk_summary':[{k:c[k] for k in ('name','version','compression','header_bytes','data_bytes','uncompressed_size')} for c in chunks]},indent=2))
