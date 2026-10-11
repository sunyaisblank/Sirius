"""Classify the actually requested code-end boundary from preserved source-bound text."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan=json.loads((W/'plan.json').read_text());trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
original=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code/query-payloads/transport.info1';join=json.loads((W/'trace-code-join.json').read_text());assert seal(original)==join['input_seals'][str(original.relative_to(ROOT))]
paths=(Path(__file__),W/'plan.json',W/'trace-code-join.json',W/'isa2-result.json',W/'rocm-gfx11-sopp-test.txt',trace,original);inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
def text(b):
 assert b[:6]==b'\x7fELF\x02\x01';position=struct.unpack_from('<Q',b,40)[0];width,n,strings=struct.unpack_from('<HHH',b,58);assert width==64
 headers=[struct.unpack_from('<IIQQQQIIQQ',b,position+i*64) for i in range(n)];s=headers[strings];names=b[s[4]:s[4]+s[5]]
 h=next(v for v in headers if names[v[0]:].split(b'\0',1)[0]==b'.text');return b[h[4]:h[4]+h[5]]
c=join['Transport']['code_object'];b=trace.read_bytes()[c['data_offset']:c['data_offset']+c['data_bytes']];assert hashlib.sha256(b).hexdigest()==c['data_sha256']
t=text(b);assert t==text(original.read_bytes()) and len(t)==237952 and hashlib.sha256(t).hexdigest()==join['Transport']['allocated_sections_and_relocations_exact']['.text']['sha256']
entry=237516;tail=t[entry:];assert t[entry-4:entry]==bytes.fromhex('0000b0bf') and tail==bytes.fromhex('00009fbf')*109
test=(W/'rocm-gfx11-sopp-test.txt').read_text();assert '# GFX11: s_code_end' in test and 'encoding: [0x00,0x00,0x9f,0xbf]' in test
r=json.loads((W/'isa2-result.json').read_text());requested={'addr':0x100209fcc,'marker_id':0};assert requested in r['first_unknown_isa_addresses'] and requested['addr']-0x1001d0000==entry
result={'source_revision':plan['source_revision'],'input_seals':inputs,'requested_unknown':requested,'source_bound_allocated_text_exact':{'bytes':len(t),'sha256':hashlib.sha256(t).hexdigest()},'entry_end_offset':entry,'tail':{'bytes':len(tail),'sha256':hashlib.sha256(tail).hexdigest(),'words':109,'word':'bf9f0000','primary_instruction':'s_code_end'},'scope':'Actual requested address is the exact entry-end in source-bound allocated text; primary GFX11 test names the identical tail word. No callback supplied this tail and no executed padding, runtime operand, full stitching, event/cost or qualification claim.'}
p=W/'code-end-byte-result.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'requested':hex(requested['addr']),'tail':result['tail'],'source_bound_text_equal':True}))
