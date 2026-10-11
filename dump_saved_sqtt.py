"""Bounded official debug-token observation of immutable saved PAL SQTT bytes."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import collections,ctypes as C,hashlib,json,os,resource,struct,threading,time
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def birth(pid):
 s=Path('/proc/'+str(pid)+'/stat').read_text();f=s[s.rfind(')')+2:].split()
 return {'pid':pid,'start_ticks':int(f[19]),'process_group':int(f[2]),'parent_pid':int(f[1])}
def write(name,v):
 p=W/name;assert not p.exists();p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
own=birth(os.getpid());parent=birth(os.getppid());start=time.monotonic()
assert sys.stdin.readline().strip()=='START'
for kind,limits in [(resource.RLIMIT_CORE,(0,0)),(resource.RLIMIT_AS,(512*1024**2,512*1024**2)),(resource.RLIMIT_CPU,(35,36)),(resource.RLIMIT_FSIZE,(16*1024**2,16*1024**2))]:resource.setrlimit(kind,limits)
stop=threading.Event()
def watch():
 while not stop.wait(.1):
  try:
   if birth(parent['pid'])!=parent or time.monotonic()-start>45:os._exit(124)
  except (FileNotFoundError,ProcessLookupError):os._exit(125)
threading.Thread(target=watch,daemon=True).start()
plan=json.loads((W/'plan.json').read_text());prep=json.loads((W/'cpu-decoder-preparation.json').read_text())
trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
entry=prep['extracted_library'];libpath=ROOT/entry['path'];assert seal(libpath)=={k:entry[k] for k in ('bytes','sha256')}
paths=[Path(__file__),W/'plan.json',W/'cpu-decoder-preparation.json',W/'decoder-preparation.json',W/'rocprof_trace_decoder.h',W/'trace_decoder_types.h',W/'spm-index-result.json',trace,libpath]
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
index=json.loads((W/'spm-index-result.json').read_text());chunks=[v for v in index['chunks'] if v['name']=='SqttData'];assert len(chunks)==1
c=chunks[0];assert c['version']==5 and c['compression']==0 and c['data_bytes']==781088
b=trace.read_bytes();raw=b[c['data_offset']:c['data_offset']+c['data_bytes']]
assert len(raw)==781088 and hashlib.sha256(raw).hexdigest()==c['data_sha256']=='cf34a1540d57146bade08c3686de5848640f5ebc35b28d6c08e4ee842800822d'
write('dump-bootstrap.json',{'source_revision':plan['source_revision'],'owner':own,'parent':parent,'input_seals':inputs,'payload':{'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'rdf_chunk_ordinal':c['ordinal']},'limits':{'wall_seconds':45,'address_space_bytes':512*1024**2,'cpu_seconds':35,'per_file_bytes':16*1024**2,'token_bytes':64*1024**2,'token_count':1000000},'format':'Each binary record is little-endian int64 time, uint32 type length, uint32 info length, then exact non-NUL type/info bytes; files preserve callback order.'})
CB=C.CFUNCTYPE(None,C.c_int64,C.c_char_p,C.c_char_p,C.c_void_p)
counts=collections.Counter();samples={};errors=[];count=0;total=0;part=0;file=None;part_bytes=0;outputs=[];first_time=None;last_time=None
def close_part():
 global file
 if file is not None:file.flush();file.close();outputs.append({'path':str(Path(file.name).relative_to(ROOT)),**seal(Path(file.name))});file=None
@CB
def debug(clock,kind,info,user):
 global count,total,part,file,part_bytes,first_time,last_time
 try:
  assert kind is not None and info is not None
  assert len(kind)<=4096 and len(info)<=16384
  count+=1;assert count<=1000000
  record=struct.pack('<qII',clock,len(kind),len(info))+kind+info
  total+=len(record);assert total<=64*1024**2
  if file is None or part_bytes+len(record)>4*1024**2:
   close_part();file=(W/('dump-tokens-%02d.bin'%part)).open('xb');part+=1;part_bytes=0
  file.write(record);part_bytes+=len(record)
  name=kind.hex();assert name in counts or len(counts)<256;counts[name]+=1
  examples=samples.setdefault(name,[])
  if len(examples)<4:examples.append({'time':clock,'info_hex':info.hex()})
  if first_time is None:first_time=clock
  last_time=clock
 except BaseException as e:
  try:
   try:errors.append({'callback':'DEBUG','error':repr(e)})
   except BaseException:pass
   try:close_part()
   except BaseException:pass
  finally:os._exit(126)
storage=C.create_string_buffer(raw,len(raw))
lib=C.CDLL(str(libpath));lib.rocprof_trace_decoder_dump_data.argtypes=[C.c_void_p,C.c_uint64,CB,C.c_void_p];lib.rocprof_trace_decoder_dump_data.restype=C.c_int
lib.rocprof_trace_decoder_get_status_string.argtypes=[C.c_int];lib.rocprof_trace_decoder_get_status_string.restype=C.c_char_p
status=lib.rocprof_trace_decoder_dump_data(C.cast(storage,C.c_void_p),len(raw),debug,None);close_part()
message=lib.rocprof_trace_decoder_get_status_string(status)
for p,v in inputs.items():assert seal(ROOT/p)==v
mapped=sorted({line.split(maxsplit=5)[5] for line in Path('/proc/self/maps').read_text().splitlines() if len(line.split(maxsplit=5))==6 and line.split(maxsplit=5)[5].startswith('/')});assert str(libpath) in mapped
result={'source_revision':plan['source_revision'],'owner':own,'input_seals':inputs,'mapped_files':mapped,'status':status,'status_string':(message or b'').decode('utf8','replace'),'callback_count':count,'type_counts_hex':dict(counts),'samples':samples,'token_files':outputs,'total_token_bytes':total,'first_callback_time':first_time,'last_callback_time':last_time,'callback_errors':errors,'elapsed_seconds':time.monotonic()-start,'scope':'Actual CPU-only public debug callback on saved bytes. Exact strings retained; undocumented grammar and event/clock mapping require separate validation. No GPU work, fabricated instrumentation, shader-cost or scientific/release verdict.'}
write('dump-result.json',result);stop.set()
print(json.dumps({'status':status,'callback_count':count,'types':len(counts),'token_bytes':total,'elapsed_seconds':result['elapsed_seconds']}),flush=True)
