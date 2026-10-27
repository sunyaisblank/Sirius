"""Check the changed raw objects before their first diagnostic execution."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R=Path.cwd();W=Path(__file__).resolve().parent
O=R/'bin/linux-gcc/native-cpu-packed-math-boundary/objects'


def pin(p):
    with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())


def inspect(p):
    text=subprocess.check_output(['objdump','-d','-M','intel','--no-show-raw-insn',str(p)],text=True)
    rows=[]
    for line in text.splitlines():
        m=re.match(r'^\s*[0-9a-f]+:\s+(\S+)(?:\s+(.*))?$',line)
        if m:
            op,args=m[1],m[2] or ''
            if op in ('{evex}','{vex}'):op,_,args=args.partition(' ')
            assert op not in ('(bad)','.byte') and not op.startswith(('vfm','vfnm'))
            rows.append((op,args))
    return rows


built=json.loads((W/'build-result.json').read_text());assert built['passed']
rows=[]
for stage,item in zip(('transport','endpoint','dopri_phase'),built['objects']):
    p=O/(stage+'.o');assert pin(p)=={k:item[k] for k in ('bytes','sha256')}
    old=R/'bin/linux-gcc/native-cpu-v4-raw-packets/objects'/(stage+'.o')
    a,b=inspect(old),inspect(p)
    additions={op for op,_ in b}-{op for op,_ in a}
    assert additions <= {'jae','or','setae','sete'},additions
    symbols=subprocess.check_output(['nm','-S','-a',str(p)],text=True)
    assert not any('RT'+op+'Store_0' in symbols for op in ('Add','Subtract','Multiply','Divide','Sqrt'))
    undefined=[line.split()[-1] for line in subprocess.check_output(['nm','-u',str(p)],text=True).splitlines()]
    allowed={'memset','memcpy','sirius_simd_fmaxf'}|{'sirius_native_v4_'+kind+'_'+op+'8' for kind in ('dense','tagged') for op in ('add','subtract','multiply','divide','sqrt')}
    assert set(undefined)<=allowed and any(n.startswith('sirius_native_v4_tagged_') for n in undefined)
    rows.append(dict(stage=stage,original=dict(path=str(old.relative_to(R)),**pin(old)),candidate=item,
        original_instructions=len(a),candidate_instructions=len(b),
        original_stack_address_operands=sum('rsp' in x or 'rbp' in x for _,x in a),
        candidate_stack_address_operands=sum('rsp' in x or 'rbp' in x for _,x in b),
        new_mnemonics=sorted(additions),undefined=undefined,store_wrappers_eliminated=True,no_FMA=True))
result=dict(passed=True,objects=rows,
    allowed_new_mnemonics='Baseline x86 jae/or/setae/sete only; all other opcodes drawn from the previously reviewed same-target raw objects.',
    explicit_compiler_private_temporary_bytes_per_store=576,
    scope='Static source/object acceptance and wrapper elimination only. Stack-address counts include LEA/frame references and are not spills or traffic. No stack-slot reuse, full-bank SSA, throughput, provider/frame or scientific/release claim.')
with (W/'code-proof.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(passed=True,raw_objects=3,no_FMA=True,store_wrappers_eliminated=True,private_temporary_bytes=576)))
