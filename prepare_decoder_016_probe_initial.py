"""Adapt the reviewed owner for the pinned 0.1.6 ABI, without invoking it."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json
W=Path(__file__).resolve().parent;ROOT=Path.cwd()
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def save(name,text):
 p=W/name;assert not p.exists();p.write_text(text)
base=W/'decode_with_code_end.py';owner=W/'run_code_end_decoder.py'
text=base.read_text().replace("W/'plan.json'","W/'decoder016-plan.json'").replace("W/'cpu-decoder-preparation.json'","W/'decoder016-preparation.json'").replace("W/'decoder-preparation.json'","W/'decoder016-elf-exports.txt'").replace("W/'rocprof_trace_decoder.h'","W/'decoder016-api.h'").replace("W/'trace_decoder_types.h'","W/'decoder016-types.h'").replace('isa-','isa016-')
text=text.replace("class Realtime(C.Structure):", "class OtherSimd(C.Structure):_fields_=[('size',C.c_uint64),('time',C.c_int64),('cycles',C.c_uint16),('wgp',C.c_uint8),('category',C.c_uint8)]\nclass Realtime(C.Structure):")
text=text.replace('assert C.sizeof(Wave)==64','assert C.sizeof(Wave)==64 and C.sizeof(OtherSimd)==24\nassert OtherSimd.time.offset==8 and OtherSimd.cycles.offset==16 and OtherSimd.wgp.offset==18 and OtherSimd.category.offset==19')
text=text.replace('isa_status_counts={}','isa_status_counts={};other=[];additional=[]')
text=text.replace('assert 0<=kind<9','assert 0<=kind<10')
needle="  return 0\n except BaseException as e:callback_error({'callback':'TRACE'"
replacement="""  elif kind in (2,6,9):
   # All three public C structures are 24 bytes in this pinned x86-64 ABI.
   length=n*24;wave_bytes+=length;assert wave_bytes<=64*1024**2
   if kind==9:
    for v in C.cast(ptr,C.POINTER(OtherSimd))[:n]:assert v.size==24
   files=[]
   for offset in range(0,length,8*1024**2):
    data=C.string_at(ptr+offset,min(8*1024**2,length-offset));p=W/('isa016-record-%02d-%04d-%03d.bin'%(kind,len(other if kind==9 else additional),offset//(8*1024**2)))
    assert not p.exists();p.write_bytes(data);files.append({'path':str(p.relative_to(ROOT)),**seal(p)})
   (other if kind==9 else additional).append({'kind':kind,'records':n,'record_bytes':24,'files':files})
  elif kind==5:assert n==0,'Undocumented DEBUG record refused'
  return 0
 except BaseException as e:callback_error({'callback':'TRACE'"""
assert text.count(needle)==1;text=text.replace(needle,replacement)
text=text.replace("'wave_array_bytes':wave_bytes", "'retained_array_bytes':wave_bytes,'other_simd_records':other,'additional_records':additional,'cu_encoding':{'shader_array_shift':7,'shader_array_mask':128,'wgp_mask':127,'scope':'Applies only to this 0.1.6 result; prior 0.1.3 CU values must not be interpreted with these macros'}")
text=text.replace('Source-bound109','Source-bound 109').replace('Actual CPU-only official decoder compatibility screen','Actual CPU-only official 0.1.6 decoder compatibility screen')
save('decode_with_decoder_016.py',text)
text=owner.read_text().replace('decode_with_code_end.py','decode_with_decoder_016.py').replace("W/'plan.json'","W/'decoder016-plan.json'").replace('isa-','isa016-')
save('run_decoder_016.py',text)
plan=json.loads((W/'plan.json').read_text());plan['purpose']='One bounded CPU-only check of official 0.1.6 on unchanged saved PAL bytes for supported other-SIMD records and timing corrections. No new GPU capture, guessed lane selection or optimization.'
plan['new_decoder']=json.loads((W/'decoder016-preparation.json').read_text())['extracted_library']
plan['input_seals']={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),base,owner,W/'plan.json',W/'decoder016-preparation.json',W/'isa-result.json',W/'isa-run.json')}
save('decoder016-plan.json',json.dumps(plan,indent=2)+'\n')
print(json.dumps({'probe':seal(W/'decode_with_decoder_016.py'),'owner':seal(W/'run_decoder_016.py'),'plan':seal(W/'decoder016-plan.json'),'native_execution':False}))
