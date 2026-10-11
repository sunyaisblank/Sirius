"""One owned, CPU-only public-ABI compatibility probe of an immutable PAL payload."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ctypes as C,hashlib,json,os,resource,threading,time
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def birth(pid):
 s=Path('/proc/'+str(pid)+'/stat').read_text();fields=s[s.rfind(')')+2:].split()
 return {'pid':pid,'start_ticks':int(fields[19]),'process_group':int(fields[2]),'parent_pid':int(fields[1])}
def write(name,v):
 p=W/name;assert not p.exists();p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
own=birth(os.getpid());parent=birth(os.getppid());start=time.monotonic()
assert sys.stdin.readline().strip()=='START'
resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2));resource.setrlimit(resource.RLIMIT_CPU,(35,36))
resource.setrlimit(resource.RLIMIT_FSIZE,(16*1024**2,16*1024**2))
stop=threading.Event()
def watch():
 while not stop.wait(.1):
  try:
   if birth(parent['pid'])!=parent or time.monotonic()-start>45:os._exit(124)
  except (FileNotFoundError,ProcessLookupError):os._exit(125)
threading.Thread(target=watch,daemon=True).start()
plan=json.loads((W/'plan.json').read_text());prep=json.loads((W/'cpu-decoder-preparation.json').read_text())
head=plan['source_revision']
trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
entry=prep['extracted_library'];libpath=ROOT/entry['path'];assert seal(libpath)=={k:entry[k] for k in ('bytes','sha256')}
paths=[Path(__file__),W/'plan.json',W/'cpu-decoder-preparation.json',W/'decoder-preparation.json',W/'rocprof_trace_decoder.h',W/'trace_decoder_types.h',W/'sdk-decode.cpp',W/'aql-threadtrace.cpp',W/'spm-index-result.json',W/'saved-isa-map.json',W/'saved-isa-summary.json',trace,libpath]
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
index=json.loads((W/'spm-index-result.json').read_text());chunks=[v for v in index['chunks'] if v['name']=='SqttData'];assert len(chunks)==1
c=chunks[0];assert c['version']==5 and c['compression']==0 and c['data_bytes']==781088
b=trace.read_bytes();raw=b[c['data_offset']:c['data_offset']+c['data_bytes']];assert hashlib.sha256(raw).hexdigest()==c['data_sha256']=='cf34a1540d57146bade08c3686de5848640f5ebc35b28d6c08e4ee842800822d'
assert len(raw)==781088
write('isa-bootstrap.json',{'source_revision':head,'owner':own,'parent':parent,'input_seals':inputs,'payload':{'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'rdf_chunk_ordinal':c['ordinal']},'limits':{'wall_seconds':45,'address_space_bytes':512*1024**2,'cpu_seconds':35,'core_bytes':0},'isa_policy':'Source-bound printed ISA, only marker zero raw VA and exact instruction boundaries; 14 archived Transport literal relocations disclosed. Missing driver/IDs rejected.'})
class PC(C.Structure):_fields_=[('addr',C.c_uint64),('marker_id',C.c_uint64)]
class Wave(C.Structure):_fields_=[('cu',C.c_uint8),('simd',C.c_uint8),('wave_id',C.c_uint8),('contexts',C.c_uint8),('r1',C.c_uint32),('r2',C.c_uint32),('r3',C.c_uint32),('begin',C.c_int64),('end',C.c_int64),('timeline_count',C.c_uint64),('instruction_count',C.c_uint64),('timeline',C.c_void_p),('instructions',C.c_void_p)]
class Realtime(C.Structure):_fields_=[('shader_clock',C.c_int64),('realtime_clock',C.c_uint64),('reserved',C.c_uint64)]
class Occupancy(C.Structure):_fields_=[('pc',PC),('time',C.c_uint64),('reserved',C.c_uint8),('cu',C.c_uint8),('simd',C.c_uint8),('wave_id',C.c_uint8),('flags',C.c_uint32)]
assert C.sizeof(Wave)==64
assert C.sizeof(PC)==16 and C.sizeof(Realtime)==24 and C.sizeof(Occupancy)==32
SE=C.CFUNCTYPE(C.c_uint64,C.POINTER(C.POINTER(C.c_uint8)),C.POINTER(C.c_uint64),C.c_void_p)
TR=C.CFUNCTYPE(C.c_int,C.c_int,C.c_void_p,C.c_uint64,C.c_void_p)
ISA=C.CFUNCTYPE(C.c_int,C.c_void_p,C.POINTER(C.c_uint64),C.POINTER(C.c_uint64),PC,C.c_void_p)
storage=(C.c_uint8*len(raw)).from_buffer_copy(raw);counts={};pairs=[];freq=[];gfx=[];infos=[];occ=[];errors=[];isa_addresses=[];isa_calls=0;se_calls=0;total=0;wave_bytes=0;waves=[];isa_status_counts={}
isa_map=json.loads((W/'saved-isa-map.json').read_text());addresses={}
for stage in isa_map['stages']:
 for v in stage['instructions']:
  addr=stage['base']+v['offset'];assert addr not in addresses;addresses[addr]=(v['instruction'].encode('ascii'),v['memory_size'])
assert len(addresses)==176335
def callback_error(v):
 try:errors.append(v)
 except BaseException:os._exit(126)
@SE
def se(buffer,size,user):
 global se_calls
 try:
  se_calls+=1
  if se_calls==1:
   assert buffer and size;buffer[0]=C.cast(storage,C.POINTER(C.c_uint8));size[0]=len(raw);return len(raw)
  assert se_calls<=3;size[0]=0;buffer[0]=C.cast(storage,C.POINTER(C.c_uint8))
  return 0
 except BaseException as e:callback_error({'callback':'SE','error':repr(e)});return 0
@ISA
def isa(instruction,memory_size,size,address,user):
 global isa_calls
 try:
  isa_calls+=1;assert isa_calls<=2000000
  value=addresses.get(address.addr) if address.marker_id==0 else None
  if value is None:
   if len(isa_addresses)<32:isa_addresses.append({'addr':address.addr,'marker_id':address.marker_id})
   code=3
  else:
   text,width=value;assert memory_size and size
   capacity=size[0];size[0]=len(text)
   if capacity<len(text):code=2
   else:
    assert instruction;C.memmove(instruction,text,len(text));memory_size[0]=width;code=0
  isa_status_counts[str(code)]=isa_status_counts.get(str(code),0)+1;return code
 except BaseException as e:callback_error({'callback':'ISA','error':repr(e)});return 1
@TR
def records(kind,ptr,n,user):
 global total,wave_bytes
 try:
  assert 0<=kind<9 and n<=1000000 and (n==0 or ptr)
  total+=n;assert total<=2000000
  counts[str(kind)]=counts.get(str(kind),0)+n
  if kind in (0,8):
   assert n<=4096
   values=list(C.cast(ptr,C.POINTER(C.c_uint64))[:n]);(gfx if kind==0 else freq).extend(values)
  elif kind==4:
   assert n<=4096;infos.extend(C.cast(ptr,C.POINTER(C.c_int))[:n])
  elif kind==7:
   assert len(pairs)+n<=50000
   for v in C.cast(ptr,C.POINTER(Realtime))[:n]:pairs.append({'shader_clock':v.shader_clock,'realtime_clock':v.realtime_clock,'reserved':v.reserved})
  elif kind==3:
   assert len(waves)+n<=4096
   for v in C.cast(ptr,C.POINTER(Wave))[:n]:
    assert v.timeline_count<=1000000 and v.instruction_count<=1000000
    parts={}
    for name,ptr,count,width in [('timeline',v.timeline,v.timeline_count,8),('instructions',v.instructions,v.instruction_count,32)]:
     length=count*width;wave_bytes+=length;assert wave_bytes<=64*1024**2 and (not length or ptr)
     files=[]
     for offset in range(0,length,8*1024**2):
      data=C.string_at(ptr+offset,min(8*1024**2,length-offset));p=W/('isa-wave-%04d-%s-%03d.bin'%(len(waves),name,offset//(8*1024**2)))
      assert not p.exists();p.write_bytes(data);files.append({'path':str(p.relative_to(ROOT)),**seal(p)})
     parts[name]={'records':count,'record_bytes':width,'files':files}
    waves.append({'cu':v.cu,'simd':v.simd,'wave_id':v.wave_id,'contexts':v.contexts,'reserved':[v.r1,v.r2,v.r3],'begin':v.begin,'end':v.end,'arrays':parts})
  elif kind==1:
   assert len(occ)+n<=100000
   for v in C.cast(ptr,C.POINTER(Occupancy))[:n]:occ.append({'addr':v.pc.addr,'marker_id':v.pc.marker_id,'time':v.time,'cu':v.cu,'simd':v.simd,'wave_id':v.wave_id,'start':bool(v.flags&1),'reserved_flags':v.flags>>1})
  return 0
 except BaseException as e:callback_error({'callback':'TRACE','kind':kind,'count':n,'error':repr(e)});return 1
lib=C.CDLL(str(libpath));lib.rocprof_trace_decoder_parse_data.argtypes=[SE,TR,ISA,C.c_void_p];lib.rocprof_trace_decoder_parse_data.restype=C.c_int
lib.rocprof_trace_decoder_get_status_string.argtypes=[C.c_int];lib.rocprof_trace_decoder_get_status_string.restype=C.c_char_p
lib.rocprof_trace_decoder_get_info_string.argtypes=[C.c_int];lib.rocprof_trace_decoder_get_info_string.restype=C.c_char_p
status=lib.rocprof_trace_decoder_parse_data(se,records,isa,None)
message=lib.rocprof_trace_decoder_get_status_string(status)
info_strings={str(i):(lib.rocprof_trace_decoder_get_info_string(i) or b'').decode('utf8','replace') for i in set(infos) if 0<=i<4}
mapped=sorted({line.split(maxsplit=5)[5] for line in Path('/proc/self/maps').read_text().splitlines() if len(line.split(maxsplit=5))==6 and line.split(maxsplit=5)[5].startswith('/')})
assert str(libpath) in mapped
for p,v in inputs.items():assert seal(ROOT/p)==v
result={'source_revision':head,'owner':own,'input_seals':inputs,'library_loaded_path':str(libpath),'mapped_files':mapped,'status':status,'status_string':(message or b'').decode('utf8','replace'),'se_callback_calls':se_calls,'isa_callback_calls':isa_calls,'isa_status_counts':isa_status_counts,'wave_records':waves,'wave_array_bytes':wave_bytes,'first_unknown_isa_addresses':isa_addresses,'record_counts':counts,'gfxip_records':gfx,'rt_frequency_hz_records':freq,'info_records':infos,'info_strings':info_strings,'realtime_records':pairs,'occupancy_records':occ,'callback_errors':errors,'elapsed_seconds':time.monotonic()-start,'scope':'Actual CPU-only official decoder compatibility screen on these saved PAL bytes. Saved source-bound assembly supplied only at proven boundaries; archived Transport literal operands include 14 unapplied loader relocations; unknown driver/IDs rejected. Nested selected-wave arrays retained exactly. No GPU work, instruction coverage, counter/event mapping, clock restoration, bottleneck or scientific/release verdict.'}
write('isa-result.json',result);stop.set()
print(json.dumps({'status':status,'status_string':result['status_string'],'records':counts,'realtime_pairs':len(pairs),'frequency':freq,'gfxip':gfx,'infos':info_strings,'isa_calls':isa_calls,'callback_errors':errors,'elapsed_seconds':result['elapsed_seconds']}),flush=True)
