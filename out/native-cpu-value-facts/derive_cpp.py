"""Reuse exact complete-value identity facts; original arithmetic is retained."""
from pathlib import Path
import hashlib,json,re
R=Path.cwd();W=Path(__file__).resolve().parent;base=R/'out/native-cpu-predicated-refinement/retained8.cpp'
original=base.read_text();source=original;changes=[]
def pin(p):return dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def function(s,name):
 m=re.search(r'__attribute__\(\(always_inline\)\) inline [^\n]+\b'+name+r'\(',s);assert m,name
 a=s.index('{',m.end());b=a+1;depth=1
 while depth:depth+=(s[b]=='{')-(s[b]=='}');b+=1
 return s[m.start():b]
def replace(a,b):
 global source
 assert source.count(a)==1,a[:100];source=source.replace(a,b,1);changes.append((a,b))
# These are precisely the original fast-domain general bodies. Extract rather
# than reconstruct any expression, and keep their operation/rounding order.
add=function(original,'fast_add_triple');mul=function(original,'fast_multiply_triple')
def general(body,name):
 start=body.index('    if(any(work)) {\n')+len('    if(any(work)) {\n')
 end=body.index('\n    }',start);text=body[start:end]
 text='\n'.join(line[4:] for line in text.splitlines())+'\n'
 assert text.count('general=finish(')==1;text=text.replace('general=finish(','return finish(',1)
 return '__attribute__((always_inline)) inline Triple '+name+'(Triple a,Triple b,U work) {\n'+text+'}\n'
new=general(add,'facts_add_general')+general(mul,'facts_multiply_general')
new+='''// Metadata bits: normalized1, known2, exact-zero4, exact-one8,
// exact-minus-one16. Complete writes set known, partial writes clear all bits.
// Legacy/unknown metadata retains only its normalization bit and recomputes
// identity facts from the same complete operand before making any decision.
__attribute__((always_inline)) inline U identity_facts(Triple a) {
    U lower=(bits(a.lo)|bits(a.tail)|bits(a.error))&u(2147483647U);
    U clean=nz(a.valid)&std::bit_cast<U>(lower==u(0));
    U h=bits(a.hi);
    return u(2)|(clean&std::bit_cast<U>((h&u(2147483647U))==u(0))&u(4))|
        (clean&std::bit_cast<U>(h==u(1065353216U))&u(8))|
        (clean&std::bit_cast<U>(h==u(3212836864U))&u(16));
}
__attribute__((always_inline)) inline U operand_facts(Triple a,U metadata) {
    U known=nz(metadata&u(2));
    if(all_on(known,u(0xffffffffU))) return metadata;
    return choose(known,metadata,(metadata&u(1))|identity_facts(a));
}
__attribute__((always_inline)) inline U negated_facts(U metadata) {
    return (metadata&~u(24))|((metadata&u(8))<<1)|((metadata&u(16))>>1);
}
'''
anchor=function(source,'register_tag')
replace(anchor,new+anchor.replace('return std::bit_cast<U>(load_words(tags+std::size_t(index)*8)==u(1));','return load_words(tags+std::size_t(index)*8);'))
replace('    U result=normalized(value)&u(1);','    U result=(normalized(value)&u(1))|identity_facts(value);')
before=function(source,'cached_add')
after='''__attribute__((always_inline)) inline Triple cached_add(Triple a,Triple b,U ta,U tb,U active) {
    U good=nz(a.valid)&nz(b.valid);
    U ma=operand_facts(a,ta),mb=operand_facts(b,tb);
    U az=nz(ma&u(4)),bz=nz(mb&u(4));
    U work=active&good&~az&~bz;Triple general=invalid();
    if(any(work)) {
        if(!all_on(nz(ta&u(1))&nz(tb&u(1)),work)) return add_triple(a,b,active);
        Triple aa=choose(work,a,constant_zero()),bb=choose(work,b,constant_zero());
        general=facts_add_general(aa,bb,work);
    }
    return choose(good,choose(az,b,choose(bz,a,general)),invalid());
}'''
replace(before,after)
before=function(source,'cached_multiply')
after='''__attribute__((always_inline)) inline Triple cached_multiply(Triple a,Triple b,U ta,U tb,U active) {
    U good=nz(a.valid)&nz(b.valid);
    U ma=operand_facts(a,ta),mb=operand_facts(b,tb);
    U is_zero=nz(ma&u(4))|nz(mb&u(4)),a_one=nz(ma&u(8)),b_one=nz(mb&u(8));
    U a_minus=nz(ma&u(16)),b_minus=nz(mb&u(16));
    U work=active&good&~is_zero&~a_one&~b_one&~a_minus&~b_minus;Triple general=invalid();
    if(any(work)) {
        U domain=nz(ta&u(1))&nz(tb&u(1))&bounded(a.hi,1518338048U)&bounded(b.hi,1518338048U);
        if(!all_on(domain,work)) return multiply_triple(a,b,active);
        Triple aa=choose(work,a,constant_zero()),bb=choose(work,b,constant_zero());
        general=facts_multiply_general(aa,bb,work);
    }
    Triple selected=choose(a_minus,negate(b),choose(b_minus,negate(a),general));
    selected=choose(a_one,b,choose(b_one,a,selected));
    selected=choose(is_zero,constant_zero(),selected);
    return choose(good,selected,invalid());
}'''
replace(before,after)
replace('    Triple result=cached_add(a,negate(b),ta,tb,active);','    Triple result=cached_add(a,negate(b),ta,negated_facts(tb),active);')
inverse=source
for a,b in reversed(changes):assert inverse.count(b)==1;inverse=inverse.replace(b,a,1)
assert inverse==original
target=W/'retained8.cpp';assert not target.exists();target.write_text(source)
proof=dict(basis=dict(path=str(base.relative_to(R)),**pin(base)),candidate=dict(path=str(target.relative_to(R)),**pin(target)),exact_inverse_restores_AK=True,
 metadata_bits=dict(normalized=1,known_identity_facts=2,exact_zero=4,exact_one=8,exact_minus_one=16),
 argument=['Each identity bit is the original RTExact numeric equality for the same five words, including signed zeros and any nonzero valid word. Known with no identity means all three exact predicates are false. Unknown metadata recomputes facts, never assumes nonidentity.',
 'Original outer valid checks and Multiply zero/a-one/b-one/a-minus/b-minus precedence remain. Add zero precedence and Subtract sign flips remain. Complete source tuples are read before any destination write. No copy is replaced by a reconstructed constant except the original Multiply zero case.',
 'Only nonidentity active valid lanes enter the original extracted fast general bodies, under the unchanged actual normalization and high2^54 Multiply bounds. Other lanes select their original identity/refusal; unknown/negative-normalization metadata retains original fallback. Inactive arithmetic operands are sanitized without changing observable words or FP controls.',
 'All original untagged arithmetic, Divide/Sqrt three corrections, final residual/radius/refusal, published five words, masks and alias handling remain. Complete tag stores now add exact facts; partial stores must invalidate all bits before reuse. No assumed postcondition normalization or literal/physical arithmetic folding.',
 'This is source reasoning plus exact extraction/inverse, not machine-checked proof, numerical/performance/frame/scientific/release acceptance. Focused original-scalar full-bank/fact/mask/alias checks and actual ISPC writers are required before a provider.' ],
 scope='One richer existing private metadata word, no extra payload/limb/register array. Cached identity tests replace repeated full-value predicates; every nonidentity expression/order is extracted from the unchanged original fast bodies. No shader/program/count/precision/FP/ISA/governor change.')
(W/'cpp-source-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(dict(source_derived=True,original_general_bodies_extracted=True,compiled=False)))
