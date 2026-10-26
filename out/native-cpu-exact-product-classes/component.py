"""Derive A/B exact product classes before the unchanged optional Double path."""
from pathlib import Path
import hashlib
import json
import sys
import zipfile

R=Path.cwd()
W=Path(__file__).resolve().parent
BASE=R/'out/native-cpu-exact-product-residual'
SEALED=R/'attestations/qualification/77c01da7/native-cpu-exact-product-residual-finite-compatibility'
PINS={
    'component.py':'11c8f32c32925343747d5bce10ed04f31028729ab54f08ad4c9a8ebf1ad24a9e',
    'retained8.cpp':'167df813bc9c56fea5c1da11aaf2c2a6120437f0f59ed31e422fa616b59a6471',
    'optional.cpp':'14ab07cafe1113506a88c70379e07548f283b8ffdf0d794595383f50049f5ab0',
    'derivation.json':'a0c3944b74e87341fa59087cc85ecbfec62f40aa1c3894cf361422de4b1e2bd4',
}
MANIFEST_SHA='ed84fbd90ec00cb221888f9e9cfc0ed48adf385ffb46acb425b86e9b02085853'
HEADER='__attribute__((always_inline)) inline void product32(F a, F b, F& hi, F& lo) {\n'
PREFIX='''    // Same RNE/gradual-underflow/DAZ-FTZ-off/masked/noncontracting scope.
    // A: <=12 significant bits per normal operand; eA+eB>=-126 keeps
    // the exact <=24-bit product normal. Both original splits are flag-free
    // with +0 low components; the original residual is +0 for every sign.
    // B: a normal power of two scales the other operand and its split low
    // exactly. eOther>=-103 keeps split components normal; eA+eB>=-103
    // keeps nonzero low*power products normal. Keep split(other) before hi
    // to retain its flags. Original residual zero has only otherLow's sign,
    // independent of the power's sign. Choose B as power when both qualify.
    // Whole-eight guards; every other vector keeps Double and fallback below.
    {
        U am=bits(a)&u(2147483647U), bm=bits(b)&u(2147483647U);
        U ae=am>>23, be=bm>>23, joint=ae+be;
        U normal=nz(ae)&nz(be)&bounded(a,1568669696U)&bounded(b,1568669696U);
        U class_a=normal&~nz((am|bm)&u(4095))&std::bit_cast<U>(joint>=u(128));
        if(all_on(class_a,u(0xffffffffU))) {
            hi=mul(a,b); lo=zero(); return;
        }
        U a_power=~nz(am&u(8388607)), b_power=~nz(bm&u(8388607));
        U other_exponent=choose(b_power,ae,be);
        U class_b=normal&(a_power|b_power)&std::bit_cast<U>(other_exponent>=u(24))&
                  std::bit_cast<U>(joint>=u(151));
        if(all_on(class_b,u(0xffffffffU))) {
            F other=choose(b_power,a,b), other_hi,other_low;
            split32(other,other_hi,other_low);
            hi=mul(a,b);
            lo=from(bits(other_low)&u(2147483648U)); return;
        }
    }
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
assert W==R/'out/native-cpu-exact-product-classes'
assert all(not (W/name).exists() for name in ('retained8.cpp','optional.cpp','derivation.json'))
manifest_path=SEALED/'MANIFEST.json'
assert pin(manifest_path)['sha256']==MANIFEST_SHA
manifest=doc(manifest_path)
assert manifest['finite_component_passed'] and not manifest['provider_adopted']
archive=R/manifest['archive']['path']
assert pin(archive)=={key:manifest['archive'][key] for key in ('bytes','sha256')}
inputs={str(manifest_path.relative_to(R)):pin(manifest_path)}
for name,sha in PINS.items():
    path=BASE/name
    assert pin(path)['sha256']==sha and pin(path)==manifest['members'][str(path.relative_to(R))]
    inputs[str(path.relative_to(R))]=pin(path)
base_proof=doc(BASE/'derivation.json')
for name,identity in base_proof['source_inputs'].items():
    assert pin(R/name)==identity,name
    inputs[name]=identity
for name,identity in base_proof['derived'].items():
    assert pin(BASE/name)==identity
member='out/native-cpu-exact-product-residual/compile-result.json'
with zipfile.ZipFile(archive) as bundle:
    compile_bytes=bundle.read(member)
assert data_pin(compile_bytes)==manifest['members'][member]
compiled=json.loads(compile_bytes)
assert compiled['exit']==0 and compiled['argv']==base_proof['compile_argv']
base_object='bin/linux-gcc/native-cpu-exact-product-residual/optional.o'
assert compiled['optional_object']==manifest['products'][base_object]
assert not compiled['original_v3_recompiled'] and not compiled['baseline_dispatcher_recompiled']
reused=[row for row in base_proof['reused_objects'] if Path(row['path']).name in ('v3.o','dispatch.o')]
assert len(reused)==2
for row in reused+[base_proof['original_v3_object']]:
    assert pin(R/row['path'])=={key:row[key] for key in ('bytes','sha256')}
compiler=Path(base_proof['original_compiler']['path'])
assert pin(compiler)=={key:base_proof['original_compiler'][key] for key in ('bytes','sha256')}
argv=list(compiled['argv'])
assert argv.count(str(BASE/'optional.cpp'))==argv.count(str(R/base_object))==1
argv[argv.index('-c')+1]=str(W/'optional.cpp')
argv[argv.index('-o')+1]=str(R/'bin/linux-gcc/native-cpu-exact-product-classes/optional.o')
assert argv[:argv.index('-c')]==compiled['argv'][:compiled['argv'].index('-c')]
assert '-march=x86-64-v4' in argv and '-ffp-contract=off' in argv

source=(BASE/'retained8.cpp').read_text()
assert source.count(HEADER)==1
candidate=source.replace(HEADER,HEADER+PREFIX,1)
assert candidate.count(HEADER+PREFIX)==1 and candidate.replace(HEADER+PREFIX,HEADER,1)==source
optional=(BASE/'optional.cpp').read_text()
old_comment='// Optional v4 only: one proved bounded product32 residual branch.'
assert optional.count(old_comment)==1
optional=optional.replace(old_comment,'// Optional v4 only: exact A/B classes, then the unchanged Double path.',1)
assert [row for row in optional.splitlines() if row.startswith('#define ')]==base_proof['unchanged_export_macro_lines']
assert '#include "retained8.cpp"' in optional
write(W/'retained8.cpp',candidate)
write(W/'optional.cpp',optional)
proof=dict(prepared_only=True,compiled=False,executed=False,source_inputs=inputs,
    derivation_script=pin(Path(__file__).resolve()),
    base_finite_manifest=dict(path=str(manifest_path.relative_to(R)),**pin(manifest_path)),
    base_compile_result=dict(archive=manifest['archive'],member=member,**data_pin(compile_bytes)),
    base_optional_object=dict(path=base_object,**compiled['optional_object']),
    derived={name:pin(W/name) for name in ('retained8.cpp','optional.cpp')},
    one_function_changed='product32',inverse_to_Double_source_exact=True,
    unchanged_Double_guard_and_body_and_original_fallback=True,
    original_v3_object=base_proof['original_v3_object'],reused_objects=reused,
    original_compiler=base_proof['original_compiler'],compile_argv=argv,
    compiler_flags_identical_to_sealed_Double=True,compile_source_and_output_paths_only_changed=True,
    unchanged_export_macro_lines=base_proof['unchanged_export_macro_lines'],
    classes=dict(
        preconditions='RNE, gradual underflow, DAZ/FTZ off, masked exceptions and no contraction; all eight lanes must satisfy each class.',
        A=dict(guard='Both operands normal, abs<=2^60, lower12 mantissa bits zero, biased exponent sum>=128.',
               proof='Each significand has at most12 bits. Both4097 splits and the normal at-most24-bit product are exact/flag-free; split lows are+0. The original final residual is+0 for all input signs.',
               result='Original float multiply high; low+0.'),
        B=dict(guard='Both operands normal and abs<=2^60; at least one has all23 fraction bits zero. Select B as power if both qualify; selected other biased exponent>=24, joint sum>=151.',
               proof='Other split lows lie on2^(eOther-23)>=2^-126. Power scaling is exact and nonzero low*power products lie on2^(eA+eB-23)>=2^-126. All omitted operations are flag-free. The original residual zero has otherLow sign only; the power sign is never XORed.',
               order='Run unchanged split32(selectedOther) before original mul(a,b), retaining exactly its sticky flags. Mixed selected sides/signs work lane-wise.',
               result='Original float multiply high; low=bitcast(bits(otherLow)&0x80000000).'),
        remainder='Literal previous Double source, including zero-residual fallback and every operation outside this primitive.'),
    scope='Source derivation only. OldDouble finite execution is provenance, not acceptance of these classes. Original AM source/v3/baseline dispatcher, five-word ABI, metadata, masks, all other functions/operations and strict flags remain unchanged. No build/workload/provider/benefit/frame/scientific/release claim; no claim of closing the~3.5x frame gap.')
write(W/'derivation.json',json.dumps(proof,indent=2)+'\n')
print(json.dumps(dict(source_only=True,one_function_changed='product32',inverse_exact=True,
    Double_and_original_fallback_unchanged=True,compiled=False,executed=False)))
