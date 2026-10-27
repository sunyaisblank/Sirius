"""Read emitted retention-multiply dependencies; do not execute the object."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import struct
import subprocess
import sys

R=Path.cwd()
W=Path(__file__).resolve().parent
B=R/'bin/linux-gcc/native-cpu-exact-product-zero-residual'
sys.dont_write_bytecode=True
sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'))
from elf_proof import Elf


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())


def listing(path):
    text=subprocess.check_output(['objdump','-d','-M','intel','--no-show-raw-insn',str(path)],text=True)
    result=[];function=None
    for line in text.splitlines():
        symbol=re.match(r'^[0-9a-f]+ <(.+)>:$',line)
        if symbol:function=symbol[1]
        match=re.match(r'^\s*([0-9a-f]+):\s+(\S+)(?:\s+(.*))?$',line)
        if not match:continue
        op,args=match[2],match[3] or ''
        if op in ('{evex}','{vex}'):op,_,args=args.partition(' ')
        assert op not in ('(bad)','.byte') and not op.startswith(('vfm','vfnm'))
        result.append(dict(address=int(match[1],16),op=op,args=args.split('#',1)[0].strip(),function=function))
    return result


def reg(operand):
    word=operand.split('{',1)[0].strip()
    vector=re.fullmatch(r'[xyz]mm(\d+)',word)
    if vector:return 'v'+vector[1]
    if re.fullmatch(r'k[0-7]',word):return word
    small={'eax':'rax','ax':'rax','al':'rax','ecx':'rcx','cx':'rcx','cl':'rcx',
           'edx':'rdx','dx':'rdx','dl':'rdx','ebx':'rbx','bx':'rbx','bl':'rbx',
           'esi':'rsi','si':'rsi','sil':'rsi','edi':'rdi','di':'rdi','dil':'rdi'}
    if word in small:return small[word]
    if re.fullmatch(r'r(?:[89]|1[0-5])d',word):return word[:-1]
    if re.fullmatch(r'r(?:ax|bx|cx|dx|si|di|[89]|1[0-5])',word):return word
    return None


path=B/'optional.o';elf=Elf(path);items=listing(path)
compiled=json.loads((W/'compile-result.json').read_text())
assert pin(path)==compiled['optional_object']
proof=json.loads((W/'derivation.json').read_text())
for row in proof['reused_objects']:
    assert pin(R/row['path'])=={key:row[key] for key in ('bytes','sha256')}
ops=Counter(x['op'] for x in items)
guards=[i for i,x in enumerate(items) if x['op']=='mov' and x['args'].endswith(',0x63ffffff')]
assert len(guards)==ops['vmulpd']==ops['vsubpd']==ops['vcvtpd2ps']==36


def is_splitter_load(index):
    item=items[index];following=items[index+1]['address']
    owner,=[s for s in elf.symbols if s['name']==item['function'] and s['section'] and s['size']]
    rows=[(location,kind,addend,target) for section,location,kind,addend,target in elf.relocations
          if section==owner['section'] and item['address']<=location<following]
    if len(rows)!=1:return False
    _,kind,addend,target=rows[0]
    if kind!=2 or not target['section']:return False
    offset=target['value']+addend+4-elf.sections[target['section']][3]
    data=elf.bytes(target['section'])
    return 0<=offset<=len(data)-4 and struct.unpack_from('<I',data,offset)[0]==0x45800800


records=[]
for guard in guards:
    first=items[guard]
    # Register scheduling may place the bound's scalar immediate on either
    # side of the multiply; select the nearest actual FP multiply, then prove
    # its splitter constant and its live dependence to the conditional.
    nearby=[i for i in range(max(0,guard-8),min(guard+10,len(items))) if items[i]['op']=='vmulps']
    assert nearby
    mul=min(nearby,key=lambda i:abs(i-guard))
    product=items[mul];assert product['function']==first['function']
    args=product['args'].split(',');assert len(args)==3
    sources={reg(s) for s in args[1:]};sources.discard(None)
    constant=None
    for i in range(mul-1,max(-1,mul-24),-1):
        item=items[i]
        if item['function']!=first['function']:break
        operand=item['args'].split(',')[0]
        target=reg(operand)
        if target in sources:
            sources.remove(target)
            if item['op']=='vbroadcastss' and is_splitter_load(i):constant=item['address'];break
    assert constant is not None,(hex(product['address']),'splitter constant not joined')
    tainted={reg(args[0])};flags=False;branch=None;trace=[]
    for item in items[mul+1:mul+23]:
        assert item['function']==first['function']
        op=item['op'];operands=item['args'].split(',')
        if op.startswith('j'):
            assert op!='jmp' and flags,(hex(item['address']),'branch not dependent on retained product')
            branch=item;break
        registers=[reg(x) for x in operands]
        if op in ('test','cmp') or op.startswith(('kortest','ktest','vptest')):
            flags=any(x in tainted for x in registers)
        elif registers and registers[0]:
            dependent=any(x in tainted for x in registers[1:])
            if op.startswith(('vpxor','vxor')) and len(registers)==3 and registers[1]==registers[2]:dependent=False
            if dependent:tainted.add(registers[0])
            else:tainted.discard(registers[0])
            if op in ('add','sub','and','or','xor','shr','shl'):flags=dependent
        trace.append(dict(address=item['address'],op=op,args=item['args']))
    assert branch is not None and first['address']<branch['address']
    records.append(dict(function=first['function'],bound_immediate=first['address'],
        splitter_broadcast=constant,splitter_bits='0x45800800',multiply=product['address'],
        conditional_branch=branch,actual_multiply_to_branch_dependency=trace))

result=dict(passed=True,optional_object=pin(path),source_derivation=pin(W/'derivation.json'),
    no_FMA=True,wide_product_sites=36,retention_multiply_sites=36,retention_dependencies=records,
    original_v3_and_baseline_dispatcher_byte_preserved=True,
    scope='Static emitted whole-object constant/operand data-dependency checks only. Each finite-bound conditional consumes the retained4097*x multiply. This is not runtime site coverage, a universal proof of the mathematical domain, a throughput measurement, or provider/frame/scientific/release acceptance.')
with (W/'code-proof.json').open('x') as stream:stream.write(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(passed=True,no_FMA=True,wide_product_sites=36,retention_multiply_dependencies=36)))
