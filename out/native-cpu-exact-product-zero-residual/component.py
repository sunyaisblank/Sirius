"""Source-only exact product/residual replacement with original zero signs/IX."""
from pathlib import Path
import hashlib
import json
import sys
import zipfile

R=Path.cwd()
W=Path(__file__).resolve().parent
OLD=R/'out/native-cpu-value-fact-stores'
FACTS=R/'out/native-cpu-exact-product-residual'
FACT_PINS={'component.py':'11c8f32c32925343747d5bce10ed04f31028729ab54f08ad4c9a8ebf1ad24a9e',
           'derivation.json':'a0c3944b74e87341fa59087cc85ecbfec62f40aa1c3894cf361422de4b1e2bd4'}
ORIGINAL='''__attribute__((always_inline)) inline void product32(F a, F b, F& hi, F& lo) {
    F ah,al,bh,bl;
    split32(a,ah,al); split32(b,bh,bl);
    F first=mul(a,b);
    hi=first;
    lo=sub(mul(al,bl),sub(sub(sub(first,mul(ah,bh)),mul(al,bh)),mul(ah,bl)));
}'''
PREFIX='''
    // RNE, DAZ/FTZ off, gradual underflow, masked exceptions, no contraction.
    // Whole-eight domain: two normals with e>=-103, |x|<=2^60, eA+eB>=-80;
    // or signed zero with such a normal; or two signed zeros. Split lows lie
    // on2^(e-23)>=2^-126; nonzero cross-products on2^(eA+eB-46)>=2^-126.
    // Thus the original TwoProduct and the binary64 residual are exact here.
    // For r=low12 bits, split low is+0 if r=0; otherwise its sign is sign(x)
    // XOR(r>2048 || r==2048 && retained bit12 odd), including exponent carry.
    // Exact residual zero has the sole nonzero split-low sign, or+0 otherwise;
    // a zero operand's sign never enters that low sign. Keep original high mul.
    // Old flags are only IX, iff some operand is noncoarse. A nonzero residual
    // makes high mul supply IX; otherwise retain one original4097*x multiply,
    // consumed by the finite guard below. Its emitted dependency must be checked.
    using D=double __attribute__((vector_size(64)));
    static_assert(sizeof(D)==64 && sizeof(double)==8 &&
                  std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits==53);
    U am=bits(a)&u(2147483647U),bm=bits(b)&u(2147483647U);
    U ae=am>>23,be=bm>>23;
    U az=~nz(am),bz=~nz(bm);
    U an=std::bit_cast<U>(ae>=u(24))&std::bit_cast<U>(am<=u(1568669696U));
    U bn=std::bit_cast<U>(be>=u(24))&std::bit_cast<U>(bm<=u(1568669696U));
    U domain=(az|an)&(bz|bn)&(az|bz|std::bit_cast<U>((ae+be)>=u(174)));
    if(all_on(domain,u(0xffffffffU))) {
        F first=mul(a,b);
        D wide_a=__builtin_convertvector(a,D),wide_b=__builtin_convertvector(b,D);
        D exact_product=wide_a*wide_b;
        D wide_high=__builtin_convertvector(first,D);
        D exact_residual=exact_product-wide_high;
        F residual=__builtin_convertvector(exact_residual,F);
        U nonzero=nz(bits(residual)&u(2147483647U));
        U ar=am&u(4095),br=bm&u(4095),al_nonzero=nz(ar),bl_nonzero=nz(br);
        bool ix_kept=true;
        if(!any(nonzero) && any(al_nonzero|bl_nonzero)) {
            F chosen=choose(al_nonzero,a,choose(bl_nonzero,b,zero()));
            F retained_split=mul(f(1166018560U),chosen);
            //4097*|chosen|<=4097*2^60<2^73; retain the actual value dependency.
            ix_kept=all_on(std::bit_cast<U>((bits(retained_split)&u(2147483647U))<u(0x64000000U)),
                           u(0xffffffffU));
        }
        if(ix_kept) {
            U a_up=std::bit_cast<U>(ar>u(2048))|
                   (std::bit_cast<U>(ar==u(2048))&nz(am&u(4096)));
            U b_up=std::bit_cast<U>(br>u(2048))|
                   (std::bit_cast<U>(br==u(2048))&nz(bm&u(4096)));
            U a_low_sign=(bits(a)^(a_up&u(2147483648U)))&al_nonzero&u(2147483648U);
            U b_low_sign=(bits(b)^(b_up&u(2147483648U)))&bl_nonzero&u(2147483648U);
            U zero_sign=(a_low_sign|b_low_sign)&(al_nonzero^bl_nonzero);
            hi=first;lo=choose(nonzero,residual,from(zero_sign));return;
        }
    }
    // Every other vector executes the literal original17-operation body below.
'''


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())


def data_pin(data):
    return dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path,text):
    with path.open('x') as stream:
        stream.write(text)


assert sys.argv[1:]==['derive']
assert W==R/'out/native-cpu-exact-product-zero-residual'
assert all(not (W/name).exists() for name in ('retained8.cpp','optional.cpp','derivation.json'))
inputs={}
for name,sha in FACT_PINS.items():
    path=FACTS/name
    assert pin(path)['sha256']==sha
    inputs[str(path.relative_to(R))]=pin(path)
facts=doc(FACTS/'derivation.json')
for name,identity in facts['source_inputs'].items():
    assert pin(R/name)==identity,name
    inputs[name]=identity
source_path=R/facts['original_source']['path']
assert source_path==OLD/'retained8.cpp'
assert pin(source_path)=={key:facts['original_source'][key] for key in ('bytes','sha256')}
source=source_path.read_text()
assert source.count(ORIGINAL)==1
replacement=ORIGINAL.replace('{\n','{'+PREFIX,1)
candidate=source.replace(ORIGINAL,replacement,1)
assert candidate.count(replacement)==1 and candidate.replace(replacement,ORIGINAL,1)==source

sealed=facts['sealed_compile_result']
archive=R/sealed['archive']['path']
assert pin(archive)=={key:sealed['archive'][key] for key in ('bytes','sha256')}
with zipfile.ZipFile(archive) as zipped:
    data=zipped.read(sealed['member'])
assert data_pin(data)=={key:sealed[key] for key in ('bytes','sha256')}
compiled=json.loads(data)
assert all(row['exit']==0 for row in compiled['commands'])
argv=list(compiled['commands'][0]['argv'])
assert '-march=x86-64-v4' in argv and '-mprefer-vector-width=256' in argv
for flag in ('-mfpmath=sse','-fexcess-precision=standard','-fno-fast-math',
             '-fno-associative-math','-fno-unsafe-math-optimizations','-frounding-math','-ffp-contract=off'):
    assert flag in argv
argv[argv.index('-c')+1]=str(W/'optional.cpp')
argv[argv.index('-o')+1]=str(R/'bin/linux-gcc/native-cpu-exact-product-zero-residual/optional.o')
reused=[row for row in facts['reused_objects'] if Path(row['path']).name in ('dispatch.o','v3.o')]
assert len(reused)==2
for row in reused+[facts['original_v3_object']]:
    assert pin(R/row['path'])=={key:row[key] for key in ('bytes','sha256')}
compiler=facts['original_compiler']
assert pin(Path(compiler['path']))=={key:compiler[key] for key in ('bytes','sha256')}

optional=(R/'out/native-cpu-v4-arithmetic/optional.cpp').read_text()
include='#include "../native-cpu-value-fact-stores/retained8.cpp"'
assert optional.count(include)==1
optional=optional.replace(include,'#include "retained8.cpp"',1)
optional=optional.replace('// Exact original math; only external names and guarded target differ.',
    '// Optional v4 only: exact residuals and original signed-zero/IX behavior.')
exports=[line for line in optional.splitlines() if line.startswith('#define ')]
assert exports==facts['unchanged_export_macro_lines'] and len(exports)==10
write(W/'retained8.cpp',candidate)
write(W/'optional.cpp',optional)
proof=dict(prepared_only=True,compiled=False,executed=False,source_inputs=inputs,
    derivation_script=pin(Path(__file__).resolve()),original_source=facts['original_source'],
    sealed_compile_result=sealed,derived={name:pin(W/name) for name in ('retained8.cpp','optional.cpp')},
    one_function_changed='product32',replacement_count=1,inverse_exact=True,
    original_fallback_body_verbatim=True,unchanged_export_macro_lines=exports,
    original_v3_object=facts['original_v3_object'],reused_objects=reused,
    original_compiler=compiler,compile_argv=argv,
    compiler_flags_identical_to_sealed_v4=True,compile_argv_source_output_paths_only_changed=True,
    proof=dict(
        preconditions='IEEE binary32/binary64; RNE, DAZ/FTZ off, gradual underflow, masked exceptions, no contraction. All eight lanes must satisfy the guard, including inactive lanes.',
        domain='Each lane: two normals with e>=-103, abs<=2^60 and eA+eB>=-80; or signed zero times one such normal; or two signed zeros. Other vectors enter the original fallback before any new FP operation.',
        exactness='Split components have lattice2^(e-23)>=2^-126. For two normals nonzero partial products have lattice2^(eA+eB-46)>=2^-126 and remain below2^123. Original TwoProduct is exact. Binary64 holds the at-most48-bit product and subtraction of the original float high exactly; binary32 residual narrowing is exact.',
        split_sign='r=absbits&4095. Original split low is+0 iff r=0; otherwise sign(input) XOR(r>2048 or r==2048 with retained bit12 odd). This includes split-high carry across a binade; both signed zeros have split low+0.',
        zero_sign='A zero residual has the sole nonzero split-low sign when exactly one low is zero; otherwise+0. For one-zero products this is the other split-low sign, independent of the zero sign. Original float multiply retains the high signed zero.',
        flags='Newly raised old flags are only IX, iff any input has nonzero r. A nonzero residual makes the retained high multiply supply IX. For all-zero residuals with a noncoarse input, one original4097*selected-noncoarse float multiply supplies IX; its absolute result is strictly<2^73 and is consumed by an actual finite guard. No artificial flag raising or FP-state writes.',
        required_code_gate='Verify the retained split multiplication and its finite-guard data dependency survive optimization. No DCE resistance is claimed from source alone. Preserve noFMA, strict flags, ISA guard and unchanged ABI/other functions.'),
    performance_scope='No hit-rate, timing or speed claim. No compiler, fixture, provider or frame execution has occurred for this component; this is not a claim to close the remaining frame gap.',
    scope='Derived directly from the exact original AO/AM retained8 source, without C2 A/B branches. Only product32 receives a prefix. Literal original fallback, every other function, external names, metadata, masks, five-word/eight-lane ABI, original v3 and baseline dispatcher remain unchanged.')
write(W/'derivation.json',json.dumps(proof,indent=2)+'\n')
print(json.dumps(dict(source_only=True,inverse_exact=True,changed_function='product32',compiled=False,executed=False)))
