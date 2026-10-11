"""Extend a proven ISA map only over source-bound Transport s_code_end padding."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan=json.loads((W/'plan.json').read_text());prior=ROOT/plan['prior']['saved-isa2-map.json']['path'];assert seal(prior)=={k:plan['prior']['saved-isa2-map.json'][k] for k in ('bytes','sha256')}
trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
original=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code/query-payloads/transport.info1';join=json.loads((W/'trace-code-join.json').read_text());assert seal(original)==join['input_seals'][str(original.relative_to(ROOT))]
paths=(Path(__file__),W/'plan.json',W/'trace-code-join.json',W/'code-end-byte-result.json',W/'rocm-gfx11-sopp-test.txt',W/'pal-code-load-header.h',W/'pal-code-load-producer.cpp',W/'pal-code-uploader.cpp',prior,trace,original);inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
def sections(b):
 assert b[:6]==b'\x7fELF\x02\x01';pos=struct.unpack_from('<Q',b,40)[0];width,n,strings=struct.unpack_from('<HHH',b,58);assert width==64
 headers=[struct.unpack_from('<IIQQQQIIQQ',b,pos+64*i) for i in range(n)];s=headers[strings];names=b[s[4]:s[4]+s[5]]
 return {names[v[0]:].split(b'\0',1)[0].decode('ascii'):{'address':v[3],'bytes':b[v[4]:v[4]+v[5]]} for v in headers}
c=join['Transport']['code_object'];b=trace.read_bytes()[c['data_offset']:c['data_offset']+c['data_bytes']];assert len(b)==c['data_bytes'] and hashlib.sha256(b).hexdigest()==c['data_sha256']
s=sections(b);os=sections(original.read_bytes());assert s['.text']==os['.text'] and s['.text']['address']==0
text=s['.text']['bytes'];assert len(text)==237952 and hashlib.sha256(text).hexdigest()==join['Transport']['allocated_sections_and_relocations_exact']['.text']['sha256']
entry=237516;tail=text[entry:];assert text[entry-4:entry]==bytes.fromhex('0000b0bf') and tail==bytes.fromhex('00009fbf')*109
test=(W/'rocm-gfx11-sopp-test.txt').read_text();assert '# GFX11: s_code_end' in test and 'encoding: [0x00,0x00,0x9f,0xbf]' in test
r=json.loads(prior.read_text());stage=next(v for v in r['stages'] if v['stage']=='Transport');assert stage['base']==0x1001d0000 and stage['entry']['size']==entry and len(stage['instructions'])==45524
assert stage['instructions'][-1]['offset']+stage['instructions'][-1]['memory_size']==entry
added=[]
for offset in range(entry,len(text),4):
 item={'offset':offset,'memory_size':4,'instruction':'s_code_end','encoded_hex':'00009fbf','provenance':'Source-bound allocated-text padding outside declared entry; primary GFX11 encoding, not claimed executed.'};stage['instructions'].append(item);added.append(stage['base']+offset)
stage['code_end_padding']={'bytes':len(tail),'sha256':hashlib.sha256(tail).hexdigest(),'instructions':len(added),'first_raw_va':added[0],'last_raw_va':added[-1]}
r['input_seals'].update(inputs);r['callback_policy']+=' Additional109 Transport allocated-text padding boundaries supply exact s_code_end; no addresses in constants/outside text or inferred auxiliary shader.'
assert sum(len(v['instructions']) for v in r['stages'])==176444
for p,v in inputs.items():assert seal(ROOT/p)==v
p=W/'saved-isa-map.json';assert not p.exists();p.write_text(json.dumps(r,separators=(',',':'))+'\n')
summary={k:v for k,v in r.items() if k!='stages'};summary['stages']=[{k:v for k,v in x.items() if k!='instructions'} for x in r['stages']];(W/'saved-isa-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
result={'source_revision':plan['source_revision'],'input_seals':inputs,'unchanged_declared_entry_records':176335,'additional_code_end_boundaries':len(added),'padding':stage['code_end_padding'],'output':{'path':str(p.relative_to(ROOT)),**seal(p)},'scope':'Exact full allocated-text/source-product match and primary GFX11 padding encoding. Future actual callback determines remaining stitching; no executed-padding, auxiliary-code, event/cost or qualification claim.'}
(W/'code-end-map-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'additional_boundaries':len(added),'requested_raw_va':hex(added[0]),'map':seal(p)}))
