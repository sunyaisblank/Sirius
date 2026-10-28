"""Gate AS against the exact sealed AR opcode vocabulary and encoding forms."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import subprocess
import sys

R=Path.cwd();W=Path(__file__).resolve().parent
O=R/'bin/linux-gcc/native-cpu-inline-private-bank/objects'
sys.dont_write_bytecode=True
sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'))
from elf_proof import Elf


def doc(p):return json.loads(p.read_text())


def pin(p):
    with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())


def checked(p,row):assert pin(p)=={k:row[k] for k in ('bytes','sha256')},p


def form(op,code):
    """Remove register-selection fields, retaining opcode/map/width/prefix form.

    This is a conservative comparison to already accepted AR encodings, not a
    general feature decoder. Any newly emitted form stops for explicit review.
    """
    i=0;prefix=[]
    while code[i] in (0xf0,0xf2,0xf3,0x2e,0x36,0x3e,0x26,0x64,0x65,0x66,0x67):
        prefix.append(code[i]);i+=1
    rex=0
    if 0x40<=code[i]<=0x4f:rex=code[i]&8;i+=1
    if code[i]==0xc5:
        return (op,'VEX',tuple(prefix),1,0,code[i+1]&7,code[i+2])
    if code[i]==0xc4:
        return (op,'VEX',tuple(prefix),code[i+1]&31,code[i+2]&128,code[i+2]&7,code[i+3])
    if code[i]==0x62:
        # Keep map, W, mandatory prefix, vector length, broadcast/rounding and
        # zero-mask form. R/X/B/vvvv and the selected k-register are operands.
        return (op,'EVEX',tuple(prefix),code[i+1]&15,code[i+2]&135,
                code[i+3]&240,bool(code[i+3]&7),code[i+4])
    opcode=[code[i]]
    if code[i]==0x0f:
        opcode.append(code[i+1])
        if code[i+1] in (0x38,0x3a):opcode.append(code[i+2])
    # These base-x86 opcode ranges encode a GP register in their low3 bits.
    if len(opcode)==1 and ((0x50<=opcode[0]<=0x5f) or (0xb0<=opcode[0]<=0xbf)):
        opcode[0]&=0xf8
    elif len(opcode)==2 and opcode[0]==0x0f and 0xc8<=opcode[1]<=0xcf:
        opcode[1]&=0xf8
    return (op,'legacy',tuple(prefix),rex,tuple(opcode))


def inspect(p):
    elf=Elf(p)
    executable={i for i,s in enumerate(elf.sections) if s[2]&4 and s[5]}
    assert executable
    section_bytes={i:elf.bytes(i) for i in executable}
    text=subprocess.check_output(['objdump','-d','-z','-M','intel','--wide','--insn-width=16',str(p)],text=True)
    rows=[];covered={i:0 for i in executable};section=None
    for line in text.splitlines():
        heading=re.match(r'^Disassembly of section (.+):$',line)
        if heading:
            matches=[i for i in executable if elf.names[i]==heading[1]]
            assert len(matches)==1,(p,heading[1]);section=matches[0];continue
        m=re.match(r'^\s*([0-9a-f]+):\s+((?:[0-9a-f]{2}\s+)+)\s*(\S.*)$',line)
        if not m:continue
        assert section in executable,(p,line)
        address=int(m[1],16);code=bytes.fromhex(m[2]);assembly=m[3]
        assert 0<len(code)<=15 and address==elf.sections[section][3]+covered[section],(p,line)
        assert section_bytes[section][covered[section]:covered[section]+len(code)]==code
        covered[section]+=len(code)
        op,_,args=assembly.partition(' ')
        while op in ('{evex}','{vex}','data16','addr32','cs','ds','ss','es','fs','gs',
                     'rep','repz','repnz','lock','notrack','rex','rex.W'):
            op,_,args=args.strip().partition(' ')
        assert op not in ('(bad)','.byte') and op
        assert not op.startswith(('vfm','vfnm','f')),(p,address,op)
        assert op not in {'ldmxcsr','stmxcsr','vldmxcsr','vstmxcsr','xsetbv','xrstor','xrstors',
                         'xsave','xsavec','xsaves','wrpkru','wrmsr','rdmsr','monitor','mwait',
                         'adcx','adox','clwb','clflushopt','rdseed','rdpid','xbegin','xend','xtest'},(p,address,op)
        rows.append(dict(section=elf.names[section],address=address,bytes=len(code),op=op,args=args,
                         encoding_form=form(op,code)))
    assert covered=={i:elf.sections[i][5] for i in executable},(p,covered)
    return elf,rows,[dict(name=elf.names[i],bytes=covered[i],sha256=hashlib.sha256(section_bytes[i]).hexdigest()) for i in sorted(executable)]


proof=doc(W/'derivation.json');built=doc(W/'build-result.json')
assert built['passed'] and built['new_raw_compilations']==3 and built['new_CXX_compilations']==0
for key in ('source_inputs','derived'):
    for name,row in proof[key].items():checked(R/name,row)
for row in proof['reference_objects']+proof['reused_objects']:checked(R/row['path'],row)
references=[];vocabulary=set();forms=set()
for row in proof['reference_objects']:
    elf,listing,sections=inspect(R/row['path'])
    references.append((elf,listing,sections))
    vocabulary|={r['op'] for r in listing}
    forms|={r['encoding_form'] for r in listing}
rows=[];removed=0
for stage,item,reference,old in zip(('transport','endpoint','dopri_phase'),built['objects'],proof['reference_objects'],references):
    p=O/(stage+'.o');checked(p,item)
    elf,new,sections=inspect(p);oldelf,original,oldsections=old
    additions={r['op'] for r in new}-vocabulary
    assert not additions,('new mnemonic requires explicit review',stage,sorted(additions))
    unexpected=[r for r in new if r['encoding_form'] not in forms]
    assert not unexpected,('new encoding form requires explicit review',stage,unexpected[:8])
    retained=[s for s in elf.symbols if 'EvaluateOneStatic_0' in s['name']]
    assert not retained,('opcode inline failed',stage,retained)
    assert not any('EvaluateOneStatic_0' in target['name'] for _,_,_,_,target in elf.relocations)
    assert not any('EvaluateOneStatic_0' in r['args'] for r in new)
    before=[s for s in oldelf.symbols if s['section'] and s['size'] and (s['info']&15)==2 and s['name'].startswith('EvaluateOneStatic_0')]
    assert len(before)=={'transport':6,'endpoint':6,'dopri_phase':4}[stage]
    removed+=len(before)
    assert not any('RT'+op+'Store_0' in s['name'] for s in elf.symbols for op in ('Add','Subtract','Multiply','Divide','Sqrt'))
    exports={s['name'] for s in elf.symbols if s['section'] and s['info']>>4 in (1,2)}
    assert len(exports)==2 and all(n.startswith('sirius_simd_v4_'+stage+'_packet') for n in exports)
    undefined={s['name'] for s in elf.symbols if not s['section'] and s['name']}
    allowed={'memset','memcpy','sirius_simd_fmaxf'}|{'sirius_native_v4_'+kind+'_'+op+'8' for kind in ('dense','tagged') for op in ('add','subtract','multiply','divide','sqrt')}
    assert undefined<=allowed and any(n.startswith('sirius_native_v4_tagged_') for n in undefined)
    shape=lambda a:dict(instructions=len(a),instruction_bytes=sum(r['bytes'] for r in a),
        stack_address_operands=sum(bool(re.search(r'\b(?:rsp|rbp)\b',r['args'])) for r in a),
        address_operand_instructions=sum('[' in r['args'] for r in a))
    rows.append(dict(stage=stage,original=reference,candidate=item,original_shape=shape(original),candidate_shape=shape(new),
        original_instructions=len(original),candidate_instructions=len(new),
        original_stack_address_operands=shape(original)['stack_address_operands'],
        candidate_stack_address_operands=shape(new)['stack_address_operands'],
        executable_sections=sections,reference_executable_sections=oldsections,
        original_mnemonics=dict(sorted(Counter(r['op'] for r in original).items())),
        candidate_mnemonics=dict(sorted(Counter(r['op'] for r in new).items())),
        new_mnemonics=sorted(additions),new_encoding_forms=[],
        removed_EvaluateOneStatic_functions=[s['name'] for s in before],EvaluateOneStatic_symbols_calls_absent=True,
        undefined=sorted(undefined),exports=sorted(exports),no_FMA=True,no_x87=True,no_FP_control_writes=True,
        all_executable_bytes_decoded=True,only_sealed_AR_instruction_forms=True))
assert removed==16
result=dict(passed=True,objects=rows,removed_EvaluateOneStatic_functions=removed,
    all_executable_sections_checked=True,unknown_or_non_v4_forms_rejected=True,
    opcode_authority=proof['parent_static_proof'],reference_objects=proof['reference_objects'],
    new_ISA_vocabulary_admitted=[],explicit_compiler_private_temporary_bytes_per_store=576,
    scope='Static changed-AS raw3 gate only. Every executable byte is decoded; mnemonic and conservative opcode/prefix/map/width forms must be observed in sealed AR raw3. No arbitrary ISA additions. All16 original opcode functions and their symbol/relocation/call references must disappear. Code/stack/address counts include frame/LEA uses and are not spills, traffic, full-bank SSA, stack reuse, timing or scientific/release acceptance.')
with (W/'code-proof.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(passed=True,raw_objects=3,removed_EvaluateOneStatic_functions=16,new_ISA_vocabulary=0)))
