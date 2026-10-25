"""Bind the original arithmetic prefix and exact cached-predicate protocol."""
from pathlib import Path
import hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent;A=R/'out/native-cpu-direct-compound'
pin=lambda p:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
a=(A/'retained8.cpp').read_text();b=(W/'retained8.cpp').read_text();assert b.startswith(a)
extra=b[len(a):]
assert extra.count('void sirius_native_tagged_')==5
assert extra.count('write_tag(tags,destination,result,active);')==5
assert 'normalized(value)&u(1)' in extra and 'if(all_on(ta&tb,work))' in extra
assert 'ta&tb&bounded(a.hi,1518338048U)&bounded(b.hi,1518338048U)' in extra
r={'passed':True,'basis':{'path':str((A/'retained8.cpp').relative_to(R)),**pin(A/'retained8.cpp')},'candidate':pin(W/'retained8.cpp'),'unchanged_AI_exact_prefix':True,
 'scope':'Five new private tagged stores. A tag is1 only for the exact unchanged AH normalized predicate on the same register words. Tags are consumed only on active lanes with both operands valid; bad operands are sanitized before the unchanged fast arithmetic and still return original canonical invalid. Multiply retains both high2^54 guards. False-negative tags use the complete original functions. Subtract uses original Add(a,Negate(b)); normalization is invariant under the three sign XORs. Divide/Sqrt use the original functions. Each active destination refreshes its tag from its complete actual result after payload stores; inactive payload/tag words remain unchanged. Both operands and tags load before any write, including aliases. No postcondition theorem, presumed input normalization, strengthened external eligibility or changed arithmetic/precision/FP controls is assumed. Independent scalar-word/std::nextafter and shared ISPC-predicate fixtures, complete writer invalidation and integration controls remain required.'}
p=W/'cpp-source-proof.json'
if p.exists():assert json.loads(p.read_text())==r
else:p.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'unchanged_AI_prefix':True,'tagged_entries':5,'compiled':False}))
