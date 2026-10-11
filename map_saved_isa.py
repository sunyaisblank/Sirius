"""Verify complete entry ISA and saved GPU allocation addresses without execution."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import collections,hashlib,json,re,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
plan=json.loads((W/'plan.json').read_text());trace=ROOT/plan['trace']['path'];assert seal(trace)=={k:plan['trace'][k] for k in ('bytes','sha256')}
index=json.loads((W/'spm-index-result.json').read_text());join=json.loads((W/'trace-code-join.json').read_text());blob=trace.read_bytes()
inputs={str(p.relative_to(ROOT)):seal(p) for p in (Path(__file__),W/'plan.json',W/'spm-index-result.json',W/'trace-code-join.json',W/'pal-code-load-header.h',W/'pal-code-load-producer.cpp',W/'pal-code-uploader.cpp',trace)}
def chunk(c,field):
 b=blob[c[field+'_offset']:c[field+'_offset']+c[field+'_bytes']];assert len(b)==c[field+'_bytes'] and hashlib.sha256(b).hexdigest()==c[field+'_sha256'];return b
load=next(c for c in index['chunks'] if c['name']=='COLoadEvent');assert load['version']==3 and load['compression']==0
n=struct.unpack('<I',chunk(load,'header'))[0];data=chunk(load,'data');assert len(data)==n*40
loads=[dict(zip(('pci_id','event_type','base','hash_lo','hash_hi','cpu_timestamp'),struct.unpack_from('<IIQQQQ',data,40*i))) for i in range(n)]
def elf_sections(b):
 assert b[:6]==b'\x7fELF\x02\x01'
 pos=struct.unpack_from('<Q',b,40)[0];width,names,strings=struct.unpack_from('<HHH',b,58);assert width==64
 headers=[struct.unpack_from('<IIQQQQIIQQ',b,pos+i*width) for i in range(names)]
 sh=headers[strings];text=b[sh[4]:sh[4]+sh[5]];sections=[]
 for i,h in enumerate(headers):
  name=text[h[0]:].split(b'\0',1)[0].decode('ascii');sections.append({'index':i,'name':name,'type':h[1],'flags':h[2],'address':h[3],'offset':h[4],'size':h[5],'link':h[6],'info':h[7],'align':h[8],'entry_size':h[9]})
 return sections
result=[]
for stage,relative,expected_count,expected_size in [('Transport','query-payloads/transport.info2',45524,237516),('Endpoint','endpoint-payloads/endpoint.info2',130811,667368)]:
 co=join[stage]['code_object'];c=next(v for v in index['chunks'] if v['ordinal']==co['ordinal']);assert c['version']==2 and c['compression']==0
 b=chunk(c,'data');assert hashlib.sha256(b).hexdigest()==co['data_sha256'];header=chunk(c,'header');pci,lo,hi=struct.unpack('<I4xQQ',header)
 matches=[v for v in loads if (v['pci_id'],v['hash_lo'],v['hash_hi'])==(pci,lo,hi)];assert len(matches)==1 and matches[0]['event_type']==0;base=matches[0]['base']
 sections=elf_sections(b);byname={v['name']:v for v in sections};text=byname['.text'];assert text['address']==0
 allocated=[];cursor=0
 for v in sections:
  if v['flags']&2:
   alignment=max(1,v['align']);assert alignment&(alignment-1)==0
   cursor=(cursor+alignment-1)&~(alignment-1);allocated.append({'name':v['name'],'gpu_offset':cursor,'size':v['size'],'alignment':alignment});cursor+=v['size']
 assert allocated[0]['name']=='.text' and allocated[0]['gpu_offset']==0
 symbols=byname['.symtab'];strings=sections[symbols['link']];s=b[strings['offset']:strings['offset']+strings['size']];entries=[]
 assert symbols['entry_size']==24
 for offset in range(symbols['offset'],symbols['offset']+symbols['size'],24):
  name,info,other,section,value,size=struct.unpack_from('<IBBHQQ',b,offset);name=s[name:].split(b'\0',1)[0].decode('ascii');entries.append({'name':name,'info':info,'section':section,'value':value,'size':size})
 entry=next(v for v in entries if v['name']=='_amdgpu_cs_main');assert entry['section']==text['index'] and entry['value']==0 and entry['size']==expected_size
 p=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code'/relative;inputs[str(p.relative_to(ROOT))]=seal(p)
 assert p.read_bytes().endswith(b'\0');encoded=bytearray();instructions=[];labels={};widths=collections.Counter();directives=[]
 for line in p.read_bytes().rstrip(b'\0').decode('ascii').splitlines():
  line=line.strip()
  if not line:continue
  if re.fullmatch(r'// _amdgpu_cs_main_pm4Crc = 0x[0-9a-f]{16}',line) or line=='_amdgpu_cs_main:' or re.fullmatch(r'\.Ltmp[0-9]+:',line) or re.fullmatch(r'\.reloc \.Ltmp[0-9]+\+(4|8), R_AMDGPU_ABS32_(HI|LO), \.rodata',line):
   directives.append({'offset':len(encoded),'text':line});continue
  m=re.fullmatch(r'_L([0-9a-fA-F]+):',line)
  if m:
   label=int(m[1],16);assert label==len(encoded) and label not in labels;labels[label]=line[:-1];continue
  asm,words=line.split(';');words=words.split();assert 1<=len(words)<=3 and all(re.fullmatch(r'[0-9a-fA-F]{8}',v) for v in words)
  offset=len(encoded);wordbytes=b''.join(struct.pack('<I',int(v,16)) for v in words);encoded.extend(wordbytes);widths[len(wordbytes)]+=1
  instructions.append({'offset':offset,'memory_size':len(wordbytes),'instruction':asm.strip(),'encoded_hex':wordbytes.hex()})
 assert len(instructions)==expected_count and len(encoded)==expected_size and bytes(encoded)==b[text['offset']:text['offset']+expected_size]
 boundaries={v['offset'] for v in instructions}
 for v in instructions:
  for label in re.findall(r'_L([0-9a-fA-F]+)',v['instruction']):assert int(label,16) in labels and int(label,16) in boundaries
 relocations=[]
 if '.rel.text' in byname:
  r=byname['.rel.text'];assert r['entry_size']==16
  for off in range(r['offset'],r['offset']+r['size'],16):
   target,info=struct.unpack_from('<QQ',b,off);symbol=entries[info>>32];typ=info&0xffffffff
   assert target in {v['offset']+4*i for v in instructions for i in range(1,v['memory_size']//4)} and typ in (1,2)
   relocations.append({'text_offset':target,'type':typ,'symbol':symbol})
 assert len(relocations)==(14 if stage=='Transport' else 0)
 result.append({'stage':stage,'base':base,'load':matches[0],'code_object':co,'entry':entry,'allocated_sections':allocated,'instruction_count':len(instructions),'encoding':{'bytes':len(encoded),'sha256':hashlib.sha256(encoded).hexdigest()},'width_counts':dict(widths),'labels':len(labels),'relocations':relocations,'directives':directives,'instructions':instructions})
for p,v in inputs.items():assert seal(ROOT/p)==v
out={'source_revision':plan['source_revision'],'input_seals':inputs,'loads':loads,'stages':result,'callback_policy':'Only marker_id zero raw VA and proven instruction boundaries. Keep archived assembly literal operands and explicitly disclose 14 Transport load relocations; no claim of live operand equality. Unknown driver address/marker IDs rejected.','scope':'Complete declared entry/source-product encoding and saved-load applicability. No native execution, event identity, clock origin, executed cost, optimization or qualification.'}
p=W/'saved-isa-map.json';assert not p.exists();p.write_text(json.dumps(out,separators=(',',':'))+'\n')
summary={k:v for k,v in out.items() if k!='stages'};summary['stages']=[{k:v for k,v in x.items() if k!='instructions'} for x in result]
(W/'saved-isa-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'stages':[(v['stage'],v['instruction_count'],v['encoding']) for v in result],'map':seal(p)}))
