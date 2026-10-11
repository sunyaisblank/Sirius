"""Adapt verified SOPP branch labels to the SDK's unsymbolized immediate syntax."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import collections,hashlib,json,re,struct
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
paths=[Path(__file__),W/'saved-isa-map.json',W/'sdk-disassembly.hpp',W/'rocm-comgr-disassembly.cpp',W/'rocm-gcn-instruction-printer.cpp',W/'rocm-amdgpu-disassembler.cpp',W/'rocm-gfx11-sopp-test.txt',W/'rocm-si-instruction-info.td']
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths};r=json.loads((W/'saved-isa-map.json').read_text());changes=[]
opcodes={'s_branch':0xbfa0,'s_cbranch_scc0':0xbfa1,'s_cbranch_scc1':0xbfa2,'s_cbranch_vccz':0xbfa3,'s_cbranch_vccnz':0xbfa4,'s_cbranch_execz':0xbfa5,'s_cbranch_execnz':0xbfa6}
for stage in r['stages']:
 count=0;negative=0;examples=[]
 for i in stage['instructions']:
  if '_L' not in i['instruction']:continue
  m=re.fullmatch(r'(s_branch|s_cbranch_(scc0|scc1|vccz|vccnz|execz|execnz)) _L([0-9a-fA-F]+)',i['instruction']);assert m and i['memory_size']==4
  word=struct.unpack('<I',bytes.fromhex(i['encoded_hex']))[0];assert word>>16==opcodes[m[1]]
  imm=word&0xffff;signed=imm if imm<32768 else imm-65536;target=i['offset']+4+4*signed;assert target==int(m[3],16)
  i['original_instruction']=i['instruction'];i['instruction']=m[1]+' '+str(imm);count+=1;negative+=signed<0
  if len(examples)<4 or (signed<0 and not any(v['signed_immediate']<0 for v in examples)):examples.append({'offset':i['offset'],'original':i['original_instruction'],'normalized':i['instruction'],'encoded_hex':i['encoded_hex'],'signed_immediate':signed,'target':target})
 assert count==(2280 if stage['stage']=='Transport' else 8980);changes.append({'stage':stage['stage'],'branches':count,'negative_immediate_branches':negative,'examples':examples})
r['input_seals'].update(inputs);r['callback_policy']+=' Only SOPP branch strings change to verified unsigned decimal imm16; encoded target equals original provider label.'
for p,v in inputs.items():assert seal(ROOT/p)==v
p=W/'saved-isa2-map.json';assert not p.exists();p.write_text(json.dumps(r,separators=(',',':'))+'\n')
summary={k:v for k,v in r.items() if k!='stages'};summary['stages']=[{k:v for k,v in x.items() if k!='instructions'} for x in r['stages']]
(W/'saved-isa2-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
result={'source_revision':r['source_revision'],'input_seals':inputs,'changes':changes,'output':{'path':str(p.relative_to(ROOT)),**seal(p)},'scope':'Byte/target-checked SDK-compatible representation only. Actual closed-decoder stitching must be observed; no new code bytes, runtime literal relocation, event/cost or scientific/release verdict.'}
(W/'branch-normalization-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'changes':[(v['stage'],v['branches'],v['negative_immediate_branches']) for v in changes],'output':seal(p)}))
