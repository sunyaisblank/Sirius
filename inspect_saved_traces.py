"""Validate the two immutable RDF indexes and their recorded counter configuration."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib,json,struct,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
assert ROOT==Path('/home/astra/.project/Sirius') and WORK==ROOT/'out/native-spm-cost'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def expected(v):return {k:v[k] for k in ('bytes','sha256')}
def git(*a):return subprocess.check_output(['git',*a],text=True).strip()
head=read(WORK/'source.json')['source_revision']
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
assert git('rev-parse','HEAD^{tree}')==git('rev-parse','b81e061^{tree}')
prior=ROOT/'attestations/native-vulkan/6be3132/native-executed-cost/manifest.json'
manifest=read(prior);assert manifest['source_revision']==head
references=read(WORK/'references.json');inputs={}
for v in references['files']:
 p=ROOT/v['path'];assert seal(p)==expected(v);inputs[v['path']]=expected(v)
inputs[str(prior.relative_to(ROOT))]=seal(prior)
for p in (Path(__file__),WORK/'source.json',WORK/'references.json'):inputs[str(p.relative_to(ROOT))]=seal(p)
rows=[]
for marker,count,sha in [('original-legacy-te.rgp',16,'ae2224a5d0b45b833d4f2ce484cd61f6488a01976d5b99959aa83749c7a9edd6'),('original-legacy-interval.rgp',20,'52420f5f1531998f48a939658165ac6c297f4d86700e68f3a283ad46458ad4e4')]:
 v=next(v for v in manifest['records'] if v['origin'].endswith('/'+marker))
 assert v['sha256']==sha and v['recovery']['kind']=='local-archive-payload'
 p=ROOT/v['recovery']['path'];assert seal(p)==expected(v);inputs[str(p.relative_to(ROOT))]=expected(v)
 b=p.read_bytes();magic,version,reserved,offset,size=struct.unpack_from('<8sIIqq',b)
 assert magic==b'AMD_RDF ' and version==3 and reserved==0 and offset>=32 and size==count*64 and offset+size==len(b)
 chunks=[];spans=[];configs=[];system=[]
 for ordinal,o in enumerate(range(offset,offset+size,64)):
  name,compression,padding,version,ho,hs,do,ds,us=struct.unpack_from('<16sB3sIqqqqq',b,o)
  assert padding==b'\0'*3 and compression==0 and us==0 and version>0
  trimmed=name.rstrip(b'\0');assert trimmed and b'\0' not in trimmed
  name=trimmed.decode('utf-8')
  assert 32<=ho<=offset and 0<=hs<=offset-ho and 32<=do<=offset and 0<=ds<=offset-do
  for start,length in ((ho,hs),(do,ds)):
   if length:spans.append((start,start+length))
  data=b[do:do+ds];header=b[ho:ho+hs]
  chunk={'ordinal':ordinal,'name':name,'version':version,'compression':compression,'header_offset':ho,'header_bytes':hs,'header_sha256':hashlib.sha256(header).hexdigest(),'data_offset':do,'data_bytes':ds,'data_sha256':hashlib.sha256(data).hexdigest()}
  if name=='TraceConfig':
   assert version==1 and hs==0 and ds==507 and data[-1:]==b'\0' and b'\0' not in data[:-1]
   config=json.loads(data[:-1].decode('utf-8'));configs.append(config);chunk['config']=config;chunk['terminal_nul_bytes']=1
  if name=='SystemInfo':
   assert version==1 and hs==0;info=json.loads(data.decode('utf-8'));system.append(info);chunk['system_info']=info
  chunks.append(chunk)
 # These concrete producer files have no overlapping or unindexed payload bytes.
 spans.sort();cursor=32
 for start,end in spans:assert start==cursor and end>start;cursor=end
 assert cursor==offset and len(configs)==len(system)==1
 config=configs[0];gpuperf=[v for v in config['sources'] if v['name']=='gpuperfexp'];assert len(gpuperf)==1
 assert gpuperf[0]['config']['spm']=={'enabled':False}
 assert gpuperf[0]['config']['sqtt']['enabled'] and gpuperf[0]['config']['sqtt']['enableInstructionTokens']
 assert not any(c['name'] in ('SpmSession','SpmCounterData','SpmDb','DerivedSpmDb') for c in chunks)
 info=system[0];assert len(info['gpus'])==1 and info['gpus'][0]['name']=='AMD Radeon 780M Graphics'
 assert info['gpus'][0]['asic']['ids']['family']==148 and info['driver']['packagingVersion'].startswith('26.10.41.01-')
 rows.append({'trace_name':marker,'protected_input':str(p.relative_to(ROOT)),**expected(next(v for v in manifest['records'] if v['origin'].endswith('/'+marker))),'index_offset':offset,'index_bytes':size,'chunks':chunks,'counter_data_available':False})
for name,value in inputs.items():assert seal(ROOT/name)==value
assert git('rev-parse','HEAD')==head==git('rev-parse','@{u}') and not git('status','--porcelain')
result={'source_revision':head,'input_seals':inputs,'traces':rows,'counter_data_available':False,'scope':'Complete concrete RDF indexes and exact TraceConfig/SystemInfo bytes inspected. Both saved traces explicitly disable SPM and lack its documented chunks; no cache/LDS/memory samples can be exported. Public format/producer references do not attest packaged binary internals or establish the exclusive reason SPM was disabled. No new scientific/device/GUI execution or optimization/render/qualification claim.'}
p=WORK/'saved-counter-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'traces':len(rows),'indexed_chunks':[len(v['chunks']) for v in rows],'counter_data_available':False,'input_versions':len(inputs)}))
