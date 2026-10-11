"""Check finite raw-instance arithmetic against stored values and pinned public formulas."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,math,re,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding='utf8-sig'))
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def digest(v):return hashlib.sha256(struct.pack('<'+str(len(v))+'d',*v)).hexdigest()
plan=read(W/'plan.json');trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
prep=read(W/'gpa-hardware-preparation.json');definitions=W/'gpa-public-counter-definitions-gfx11.txt'
assert seal(definitions)=={'bytes':51492,'sha256':'6d7cb2f0c542102d30ece3259d7caa6ea2c0a492339342d9b9f78c1f3cc20419'}
paths=[Path(__file__),W/'plan.json',trace,W/'spm-index-result.json',W/'raw-spm-result.json',W/'derived-spm-result.json',W/'gpa-hardware-preparation.json',definitions]+[ROOT/v['path'] for v in prep['files']]
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
for v in prep['files']:assert seal(ROOT/v['path'])=={k:v[k] for k in ('bytes','sha256')}
index=read(W/'spm-index-result.json');raw=read(W/'raw-spm-result.json');meta=read(W/'derived-spm-result.json');b=trace.read_bytes();n=raw['session']['timestamps'];assert n==11911
def data(c):
 assert c['compression']==0
 d=b[c['data_offset']:c['data_offset']+c['data_bytes']];assert len(d)==c['data_bytes'] and hashlib.sha256(d).hexdigest()==c['data_sha256'];return d
streams={}
for c in index['chunks']:
 if c['name']=='SpmCounterData':
  pci,block,instance,event,width=struct.unpack_from('<5I',b,c['header_offset']);assert width==2 and pci==raw['session']['pci_id']
  key=(block,event);assert instance not in streams.setdefault(key,{})
  streams[key][instance]=struct.unpack('<'+str(n)+'H',data(c))
assert len(streams)==27 and sum(map(len,streams.values()))==183
sums={k:[float(sum(v[i] for v in instances.values())) for i in range(n)] for k,instances in streams.items()}
maximum={k:[float(max(v[i] for v in instances.values())) for i in range(n)] for k,instances in streams.items()}
stored={}
for v in meta['counters']:stored[v['counter_index']]=struct.unpack_from('<'+str(n)+'d',data(index['chunks'][v['chunk_ordinal']]))
assert set(stored)==set(range(33))
public={}
for part in definitions.read_text().split('name=')[1:]:
 lines=part.splitlines();equations=[v[4:] for v in lines if v.startswith('eqn=')]
 if len(equations)==1:public[lines[0]]=equations[0]
names={}
for f in prep['files']:
 group=Path(f['path']).stem.removeprefix('gpa-hardware-').upper()
 for num,name,prefix in re.findall(r'GpaHardwareCounterDesc\{(\d+), "([^"]+)", "([^"]+)"',(ROOT/f['path']).read_text()):
  if prefix==group+('' if group=='CPF' else '0'):names[(group,int(num))]=name
binding={0:'CPF',8:'TA',10:'TCP',29:'GCEA',35:'GL1C',38:'GL2C',49:'SQWGP'}
numeric_bindings=[]
for block,event in sorted(sums):
 if block in binding:
  assert (binding[block],event) in names
  numeric_bindings.append({'raw_block':block,'event':event,'public_instance_zero_name':names[(binding[block],event)],'observed_instances':sorted(streams[(block,event)])})
def add(*v):return [sum(a) for a in zip(*v)]
def subtract(a,z):return [x-y for x,y in zip(a,z)]
def scale(a,k):return [x*k for x in a]
def percentage(a,z):return [(x/y)*100 if y else 0.0 for x,y in zip(a,z)]
def cache_percentage(miss,requests):return [(1-x/y)*100 if y else 0.0 for x,y in zip(miss,requests)]
expected={};refs={};recipes={}
def record(idx,name,values,recipe):
 assert name in public and len(values)==n;expected[idx]=values;refs[idx]=name;recipes[idx]=recipe
for base,prefix,events in [(0,'Inst',(270,271,272)),(4,'Scalar',(294,295,296))]:
 hits=sums[49,events[0]];miss=add(sums[49,events[1]],sums[49,events[2]]);requests=add(hits,sums[49,events[1]],sums[49,events[2]])
 for offset,name,values,recipe in [(0,prefix+'CacheHit',percentage(hits,requests),'hits / (hits + misses + duplicate misses) * 100; zero denominator -> 0'),(1,prefix+'CacheRequestCount',requests,'sum hits + sum misses + sum duplicate misses'),(2,prefix+'CacheHitCount',hits,'sum hits'),(3,prefix+'CacheMissCount',miss,'sum misses + sum duplicate misses')]:record(base+offset,name,values,recipe)
for pct,request_idx,hit_idx,miss_idx,prefix,key1,key2 in [(8,9,10,13,'L0',(10,9),(10,17)),(14,15,16,17,'L1',(35,14),(35,18)),(18,19,20,21,'L2',(38,3),(38,43))]:
 requests=sums[key1];miss=sums[key2]
 record(pct,prefix+'CacheHit',cache_percentage(miss,requests),'(1 - sum misses / sum requests) * 100; zero denominator -> 0')
 record(request_idx,prefix+'CacheRequestCount',requests,'sum requests')
 record(hit_idx,prefix+'CacheHitCount',subtract(requests,miss),'sum requests - sum misses')
 record(miss_idx,prefix+'CacheMissCount',miss,'sum misses')
gpu=sums[0,24];busy=maximum[8,15];stalled=maximum[10,39]
record(24,'GPUBusyCycles',gpu,'raw CPF busy cycles instance 0')
record(28,'MemUnitBusyCycles',busy,'maximum TA busy cycles over the 12 observed instances')
record(29,'MemUnitStalledCycles',stalled,'maximum TCP request-stall cycles over the 12 observed instances')
record(26,'MemUnitBusy',percentage(busy,gpu),'maximum TA busy / CPF busy * 100; zero denominator -> 0')
record(27,'MemUnitStalled',percentage(stalled,gpu),'maximum TCP request stall / CPF busy * 100; zero denominator -> 0')
lds=[x/24 for x in sums[49,256]]
record(25,'CSLDSBankConflictCycles',lds,'sum observed SQWGP bank-conflict cycles / 24; divisor fitted to this recorded stream, not independently bound NUM_SIMDS metadata')
record(23,'CSLDSBankConflict',percentage(lds,gpu),'sum observed SQWGP bank-conflict cycles / 24 / CPF busy * 100; zero denominator -> 0; same fitted-divisor limitation')
fetch=add(*(scale(sums[38,event],size) for event,size in ((99,32),(100,64),(101,96),(102,128))))
record(30,'FetchSize',fetch,'sum 32/64/96/128-byte GL2C read requests multiplied by their sizes')
write=add(scale(subtract(sums[38,83],sums[38,85]),32),scale(sums[38,85],64))
record(12,'WriteSize',write,'(sum MC_WRREQ - sum EA_WRREQ_64B) * 32 + sum EA_WRREQ_64B * 64')
record(31,'LocalVidMemBytes',scale(sums[29,55],32),'sum observed GCEA DRAM size requests * 32')
checks=[]
for idx in sorted(expected):
 actual=stored[idx];values=expected[idx];diff=[i for i,(x,y) in enumerate(zip(values,actual)) if struct.pack('<d',x)!=struct.pack('<d',y)]
 checks.append({'derived_index':idx,'display_name':meta['counters'][idx]['display_name'],'public_counter_name':refs[idx],'public_expression':public[refs[idx]],'observed_recipe':recipes[idx],'samples':n,'binary64_different_samples':len(diff),'first_different_samples':diff[:8],'maximum_absolute_error':max(abs(x-y) for x,y in zip(values,actual)),'expected_binary64_sha256':digest(values),'stored_binary64_sha256':digest(actual),'minimum':min(actual),'maximum':max(actual),'nonzero_samples':sum(v!=0 for v in actual)})
assert len(checks)==30
for p,v in inputs.items():assert seal(ROOT/p)==v
result={'source_revision':plan['source_revision'],'input_seals':inputs,'trace':seal(trace),'public_gpa_revision':prep['revision'],'numeric_name_bindings':numeric_bindings,'counter_checks':checks,'unreconstructed_indices':[11,22,32],'limits':['Finite all-sample arithmetic agreement and pinned public formulas only; no packaged GPA binary build/formula/filter attestation or complete hardware inventory.','LDS divisor 24 is an observed fitting parameter; NUM_SIMDS applicability remains separately unbound.','Ray zero streams and writer-placeholder PCIe zero stream not treated as unique formula reconstruction.','No SPM-to-SQTT/RGP event origin mapping, exclusive ownership, ordinary performance benefit or scientific/release qualification.']}
p=W/'public-formula-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'checked_streams':len(checks),'samples_per_stream':n,'bit_exact_streams':sum(v['binary64_different_samples']==0 for v in checks),'differences':[v for v in checks if v['binary64_different_samples']]}))
