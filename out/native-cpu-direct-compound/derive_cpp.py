"""Verify the unchanged AH prefix and exact original subtraction composition."""
from pathlib import Path
import hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;A=R/'out/native-cpu-normalized-arithmetic'
pin=lambda p:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
a=(A/'retained8.cpp').read_text();b=(W/'retained8.cpp').read_text();assert b.startswith(a)
addition=b[len(a):];assert addition.count('void sirius_native_dense_subtract8(')==1
assert 'Triple a=read(scratch,a_index), b=read(scratch,b_index);' in addition
assert 'write(scratch,destination,add_triple(a,negate(b),active),active);' in addition
original=R/'bin/linux-gcc/native-cpu-static-program/typed-batch/transport.native.hpp'
s=original.read_text();start=s.index('RetainedTriple_0 RTSubtract_0(');start=s.index('{',start);depth=1;end=start+1
while depth:depth+=(s[end]=='{')-(s[end]=='}');end+=1
body=s[start:end];assert 'RTNegate_0' in body and 'RTAdd_0' in body
result={'passed':True,'basis':{'path':str((A/'retained8.cpp').relative_to(R)),**pin(A/'retained8.cpp')},'candidate':pin(W/'retained8.cpp'),'unchanged_AH_exact_prefix':True,'original_scalar':{'path':str(original.relative_to(R)),**pin(original)},'original_subtract_body':body,'scope':'Only one new direct-bank entry composes unchanged AH read/negate/add/write. Both operands snapshot before writes; masks, inactive words and alias semantics retained. Exact arithmetic semantics use original scalar subtraction composition, not a new formula. All prior AH source byte-identical. Finite fixture and actual linked-code checks still required.'}
p=W/'cpp-source-proof.json'
if p.exists():assert json.loads(p.read_text())==result
else:p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'unchanged_AH_prefix':True,'new_subtract_entry':1,'compiled':False}))
