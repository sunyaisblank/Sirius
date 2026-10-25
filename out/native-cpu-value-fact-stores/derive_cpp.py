"""Use complete cached value facts for original identity-only bank stores."""
from pathlib import Path
import hashlib,json,re
R=Path.cwd();W=Path(__file__).resolve().parent;base=R/'out/native-cpu-value-facts/retained8.cpp';original=base.read_text();source=original;changes=[]
def pin(p):return dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def replace(a,b):
 global source
 assert source.count(a)==1,a[:100];source=source.replace(a,b,1);changes.append((a,b))
replace('    return u(2)|(clean&','    return u(2)|(nz(a.valid)&u(32))|(clean&')
replace('// exact-minus-one16. Complete writes set known, partial writes clear all bits.','// exact-minus-one16, nonzero-valid32. Known certifies every fact, including\n// the exact normalization bit. Complete writes set known; partial writes clear all bits.')
anchor='} // namespace\n\nextern "C" __attribute__((visibility("hidden")))'
# The original anonymous namespace ends earlier; insert into the final cached
# namespace immediately after cached_multiply, before tagged exports.
marker='\n}\n\nextern "C" __attribute__((visibility("hidden")))\nvoid sirius_native_tagged_add8'
assert source.count(marker)==1
helpers='''
__attribute__((always_inline)) inline void store_known_facts(std::uint32_t* tags,
    std::uint32_t index,U metadata,U active) {
    auto* p=tags+std::size_t(index)*8;
    if(_mm256_movemask_ps(std::bit_cast<__m256>(active))!=255)
        metadata=choose(active,metadata,load_words(p));
    std::memcpy(p,&metadata,sizeof(metadata));
}
// Only complete, coherent metadata enters this path. Every active lane must
// select an original identity/refusal; otherwise no store occurs here.
__attribute__((always_inline)) inline bool store_identity(std::uint32_t* scratch,
    std::uint32_t* tags,std::uint32_t destination,std::uint32_t ai,std::uint32_t bi,
    U ma,U mb,U active,unsigned operation) {
    if(!any(active)) return true;
    if(!all_on(nz(ma&u(2))&nz(mb&u(2)),active)) return false;
    U live=active&nz(ma&u(32))&nz(mb&u(32));
    U az=nz(ma&u(4)),bz=nz(mb&u(4));
    U zero_mask=u(0),copy_a=u(0),copy_b=u(0),negative_a=u(0),negative_b=u(0);
    if(operation==4) {
        zero_mask=live&(az|bz);live=live&~zero_mask;
        copy_b=live&nz(ma&u(8));live=live&~copy_b;
        copy_a=live&nz(mb&u(8));live=live&~copy_a;
        negative_b=live&nz(ma&u(16));live=live&~negative_b;
        negative_a=live&nz(mb&u(16));live=live&~negative_a;
    } else {
        copy_b=live&az;live=live&~copy_b;
        copy_a=live&bz;live=live&~copy_a;
        if(operation==3) {negative_b=copy_b;copy_b=u(0);}
    }
    if(any(live)) return false;
    Triple a=invalid(),b=invalid();
    // Snapshot every selected operand before writing, including destination
    // aliases. Zero/refusal-only packets load no payload operands.
    if(any(copy_a|negative_a)) a=read(scratch,ai);
    if(any(copy_b|negative_b)) b=read(scratch,bi);
    Triple result=choose(zero_mask,constant_zero(),choose(copy_b,b,choose(copy_a,a,
        choose(negative_b,negate(b),choose(negative_a,negate(a),invalid())))));
    U metadata=choose(zero_mask,u(39),choose(copy_b,mb,choose(copy_a,ma,
        choose(negative_b,negated_facts(mb),choose(negative_a,negated_facts(ma),u(2))))));
    write(scratch,destination,result,active);
    store_known_facts(tags,destination,metadata,active);
    return true;
}
'''
replace(marker,helpers+marker)
for operation,name in [(2,'add'),(4,'multiply'),(3,'subtract')]:
 start=source.index('void sirius_native_tagged_'+name+'8(');end=source.index('\n}',start)+2;before=source[start:end]
 a='    Triple a=read(scratch,a_index), b=read(scratch,b_index);\n    U active=active_mask(mask),ta=register_tag(tags,a_index),tb=register_tag(tags,b_index);'
 b='    U active=active_mask(mask),ta=register_tag(tags,a_index),tb=register_tag(tags,b_index);\n    if(store_identity(scratch,tags,destination,a_index,b_index,ta,tb,active,'+str(operation)+')) return;\n    Triple a=read(scratch,a_index), b=read(scratch,b_index);'
 assert before.count(a)==1;replace(before,before.replace(a,b,1))
inverse=source
for a,b in reversed(changes):assert inverse.count(b)==1;inverse=inverse.replace(b,a,1)
assert inverse==original
p=W/'retained8.cpp';assert not p.exists();p.write_text(source)
proof=dict(basis=dict(path=str(base.relative_to(R)),**pin(base)),candidate=dict(path=str(p.relative_to(R)),**pin(p)),exact_inverse_restores_basic_value_facts=True,
 metadata_bits=dict(normalized=1,complete_current_facts=2,exact_zero=4,exact_one=8,exact_minus_one=16,nonzero_valid=32),
 argument=['Complete metadata now certifies all facts and the exact normalization bit. Unknown/legacy/false-negative metadata cannot enter the identity-store path. Partial writes must clear all bits. The valid bit means original valid!=0, not canonical valid1.',
 'Only packets whose every active lane selects an original identity or invalid-input refusal return early. Multiply preserves zero/a-one/b-one/a-minus/b-minus precedence. Add preserves a-zero before b-zero; Subtract uses the same Add precedence with bitwise negated b. Source valid words and zero signs are copied exactly when the original copies them.',
 'All selected source payloads and metadata are snapshotted before destination stores. Zero/refusal packets require no payload loads. Copy/negate metadata equals freshly computed complete facts: negation preserves normalization/valid/zero and swaps exact-one/minus-one. Canonical zero facts39 and canonical invalid facts2 are exact.',
 'Any active general lane, or any active unknown metadata, falls through without a store to the previously checked value-facts implementation. All original untagged arithmetic and nonidentity expressions, Divide/Sqrt refinements, masks/public words/order remain. Fewer unused operations may change unobserved intermediate FP exception flags; original controls and worker full-state restoration remain.',
 'Source reasoning/exact inverse only, not machine-checked proof or performance/science/release acceptance. Complete original-scalar bank/fact/mask/alias checks and actual new writer regressions precede a provider.' ],
 scope='Existing private metadata word, no extra published limb or context array. Original identity/refusal results can be stored directly from coherent facts. No formula, precision eligibility, FP/ISA, canonical9380tuples, fields/RHS/trials/quartics, requests/count/cache/governor change.')
(W/'cpp-source-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(dict(identity_bank_stores_derived=True,compiled=False)))
