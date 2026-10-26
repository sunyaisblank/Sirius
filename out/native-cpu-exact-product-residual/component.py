"""Derive only one bounded exact-residual branch in the optional-v4 source."""
from pathlib import Path
import hashlib
import json
import sys
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
OLD = R / 'out/native-cpu-value-fact-stores'
V4 = R / 'out/native-cpu-v4-arithmetic'
SEALED = R / 'attestations/qualification/fecafd04/native-cpu-v4-arithmetic-finite-compatibility'
PINS = {
    'out/native-cpu-value-fact-stores/retained8.cpp': '60455d362b1c98b30e9975c35011c6a193ffcf98a5872ddb45036c0e5364961d',
    'out/native-cpu-value-fact-stores/prepared-cpp.json': '0631a88a4f7deb207073ee2c4219526ef77f2b668e128d6587cae24226f5c053',
    'out/native-cpu-value-fact-stores/cpp-build-result.json': '60c215d747ce87f5d2127f03d2ecb8937429f11c0a611472e7fa2aa30719945d',
    'out/native-cpu-v4-arithmetic/optional.cpp': 'bcd7804fc7935de78846bf6317d9545b7df78afa934957d74691e74927b028d5',
    'out/native-cpu-v4-arithmetic/dispatch.cpp': 'ab1512b576709665846a60b5e68f0c4dc9dbe86c08a84e794c18622037117445',
    'out/native-cpu-v4-arithmetic/component.py': '73ec96434f84b870382c5dbbf8d93efb39d80f2963ec08d62f41cacee7d781f0',
    'out/native-cpu-v4-arithmetic/derivation.json': '4bfaf91e35be759cdec724798f6a08db4e1438fae7e355b69b527d3ab63bb803',
    'attestations/qualification/fecafd04/native-cpu-v4-arithmetic-finite-compatibility/MANIFEST.json': '70220651fa964503111d9ada142217e12b884da836b1568d9a9b9ac017df5e19',
    'src/sirius/kernels/retained_pair.slang': '395dd782f3309b22f0fb481a2bf2a898489e700043c248c0fad2e53dd2fbe304',
    'out/native-cpu-slang-bridge/shaders/retained_scalar.slang': '5e74f8059512b417c2236f360264cd6f980603e06ae5e34bed5b39cabefe574e',
}

ORIGINAL = '''__attribute__((always_inline)) inline void product32(F a, F b, F& hi, F& lo) {
    F ah,al,bh,bl;
    split32(a,ah,al); split32(b,bh,bl);
    F first=mul(a,b);
    hi=first;
    lo=sub(mul(al,bl),sub(sub(sub(first,mul(ah,bh)),mul(al,bh)),mul(ah,bl)));
}'''

PREFIX = '''
    // RNE, gradual underflow, DAZ/FTZ off, masked exceptions, no contraction.
    // For normal |a|,|b|<=2^60, eA,eB>=-103 and eA+eB>=-80, the original 4097
    // split is exact and every nonzero split component is normal. Every nonzero
    // subsequent product lies on the 2^(eA+eB-46) lattice, hence is normal;
    // no splitting/product overflow is possible. TwoProduct therefore returns
    // hi=RN32(a*b), lo=a*b-hi exactly (Ogita/Rump/Oishi 2005, Alg3.3/Thm3.4).
    // Binary64 holds the 48-bit product and its residual exactly; narrowing a
    // nonzero residual is exact. Keep the original float high multiplication.
    // All-eight nonzero residuals also retain its inexact sticky flag. Smaller
    // normal inputs can have exact subnormal split lows/denormal-operand flags;
    // they retain the original body. Zero residuals use the original body:
    // 0x3f800fff*2 can have original lo=-0.
    // No active-mask/ABI/caller/other primitive or rounding policy changes.
    using D = double __attribute__((vector_size(64)));
    static_assert(sizeof(D)==64 && sizeof(double)==8 &&
                  std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits==53);
    U am=bits(a)&u(2147483647U), bm=bits(b)&u(2147483647U);
    U ae=am>>23, be=bm>>23;
    U safe=std::bit_cast<U>(ae>=u(24))&std::bit_cast<U>(be>=u(24))&
           bounded(a,1568669696U)&bounded(b,1568669696U)&
           std::bit_cast<U>((ae+be)>=u(174)); // (eA+127)+(eB+127)>=174.
    if(all_on(safe,u(0xffffffffU))) {
        F first=mul(a,b);
        D wide_a=__builtin_convertvector(a,D), wide_b=__builtin_convertvector(b,D);
        D exact_product=wide_a*wide_b;
        D wide_high=__builtin_convertvector(first,D);
        D exact_residual=exact_product-wide_high;
        F residual=__builtin_convertvector(exact_residual,F);
        if(all_on(nz(bits(residual)&u(2147483647U)),u(0xffffffffU))) {
            hi=first; lo=residual; return;
        }
    }
    // Every other vector executes the unchanged original chain below.
'''


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def byte_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def doc(path):
    return json.loads(path.read_text())


def named(row):
    path = R / row['path']
    assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}, path
    return path


def derive():
    assert W == R / 'out/native-cpu-exact-product-residual'
    assert all(not (W / name).exists() for name in ('retained8.cpp', 'optional.cpp', 'derivation.json'))
    basis = {}
    for name, expected in PINS.items():
        identity = pin(R / name)
        assert identity['sha256'] == expected, name
        basis[name] = identity
    manifest = doc(SEALED / 'MANIFEST.json')
    archive = named(manifest['archive'])
    member = 'out/native-cpu-v4-arithmetic/compile-result.json'
    with zipfile.ZipFile(archive) as source:
        data = source.read(member)
    assert byte_pin(data) == manifest['members'][member]
    compiled = json.loads(data)
    assert all(row['exit'] == 0 for row in compiled['commands'])
    assert not compiled['driver_recompiled'] and not compiled['original_v3_recompiled']
    reused = []
    for name in ('v3.o', 'dispatch.o', 'optional.o'):
        path = R / 'bin/linux-gcc/native-cpu-v4-arithmetic' / name
        relative = str(path.relative_to(R))
        assert pin(path) == manifest['products'][relative] == compiled['objects'][path.stem]
        reused.append(dict(path=relative, **pin(path)))
    original_object = compiled['original_v3']
    named(original_object)
    prepared = doc(OLD / 'prepared-cpp.json')
    built = doc(OLD / 'cpp-build-result.json')
    assert built['passed']
    old_compile, = [row for row in prepared['compilations'] if row['name'] == 'retained8']
    assert next(row['object'] for row in built['compilations'] if row['name'] == 'retained8') == original_object
    command = list(compiled['commands'][0]['argv'])
    assert command.count('-march=x86-64-v4') == command.count('-mprefer-vector-width=256') == 1
    for flag in ('-mfpmath=sse', '-fexcess-precision=standard', '-fno-fast-math',
                 '-fno-associative-math', '-fno-unsafe-math-optimizations', '-frounding-math', '-ffp-contract=off'):
        assert flag in command, flag
    flags = command[:command.index('-c')]
    original_flags = old_compile['argv'][:old_compile['argv'].index('-c')]
    expected_flags = list(original_flags)
    expected_flags[expected_flags.index('-march=x86-64-v3')] = '-march=x86-64-v4'
    expected_flags.append('-mprefer-vector-width=256')
    assert flags == expected_flags
    compiler = Path(command[0])
    assert pin(compiler) == prepared['inputs'][str(compiler)]
    command[command.index('-c') + 1] = str(W / 'optional.cpp')
    command[command.index('-o') + 1] = str(R / 'bin/linux-gcc/native-cpu-exact-product-residual/optional.o')

    original = (OLD / 'retained8.cpp').read_text()
    assert original.count(ORIGINAL) == 1
    replacement = ORIGINAL.replace('{\n', '{' + PREFIX, 1)
    candidate = original.replace(ORIGINAL, replacement, 1)
    assert candidate.count(replacement) == 1 and candidate.replace(replacement, ORIGINAL, 1) == original
    optional = (V4 / 'optional.cpp').read_text()
    include = '#include "../native-cpu-value-fact-stores/retained8.cpp"'
    assert optional.count(include) == 1
    optional = optional.replace(include, '#include "retained8.cpp"', 1)
    optional = optional.replace('// Exact original math; only external names and guarded target differ.',
                                '// Optional v4 only: one proved bounded product32 residual branch.')
    exports = [line for line in optional.splitlines() if line.startswith('#define ')]
    assert len(exports) == 10 and exports == [line for line in (V4 / 'optional.cpp').read_text().splitlines()
                                             if line.startswith('#define ')]
    for name, text in (('retained8.cpp', candidate), ('optional.cpp', optional)):
        with (W / name).open('x') as stream:
            stream.write(text)
    proof = dict(
        prepared_only=True, compiled=False, executed=False,
        source_inputs=basis, derivation_script=pin(Path(__file__).resolve()),
        stable_component_manifest=dict(path=str((SEALED / 'MANIFEST.json').relative_to(R)),
                                       **pin(SEALED / 'MANIFEST.json')),
        sealed_compile_result=dict(archive=manifest['archive'], member=member, **byte_pin(data)),
        original_source=dict(path=str((OLD / 'retained8.cpp').relative_to(R)), **pin(OLD / 'retained8.cpp')),
        derived={name: pin(W / name) for name in ('retained8.cpp', 'optional.cpp')},
        one_function_changed='product32', replacement_count=1, inverse_exact=True,
        original_fallback_body_verbatim=True, unchanged_export_macro_lines=exports,
        original_v3_object=original_object, reused_objects=reused,
        original_compiler=dict(path=str(compiler), **pin(compiler)),
        compile_argv=command, compiler_flags_identical_to_sealed_v4=True,
        compile_argv_source_output_paths_only_changed=True,
        proof=dict(
            preconditions='IEEE binary32/binary64, RNE, DAZ/FTZ off, masked exceptions, strict noncontracting vector arithmetic.',
            normal_input_domain='All eight lanes: unbiased eA,eB in[-103,60], |a|,|b|<=2^60, eA+eB>=-80. Individual biased exponents>=24 and their sum>=174.',
            split='Original4097 split is exact with no scaling overflow/underflow. Each input is a multiple of2^(e-23)>=2^-126, as are its split components; all nonzero components are normal. Smaller normal inputs can have exact subnormal split lows (for example0x00800001); the stricter individual exponent guard retains the original path and its denormal-operand flag.',
            partial_products='Each split component is a multiple of its original2^(e-23) lattice. Each nonzero cross-product is a multiple of2^(eA+eB-46)>=2^-126. All five original product magnitudes remain below2^123, so none overflows or rounds through underflow.',
            equality='Original TwoProduct has high=RN32(a*b), low=a*b-high exactly in this domain. A binary64 significand holds the exact at-most48-bit binary32 product. Its subtraction of the same original float high is exact; the resulting nonzero residual is exactly representable as binary32. Thus both original primitive words are identical.',
            zero='Any zero residual in any lane takes the entire original chain, retaining all zero signs. Example0x3f800fff*2 has original low=-0 while wide exact cancellation gives+0.',
            exceptions='Every fast lane has nonzero residual, so the retained original high multiplication already raises inexact. The skipped original chain can raise no other exception in this domain; all wide operations/conversions are exact. Any fallback executes the unchanged original chain.',
            inactive_and_unknown='Whole-eight guard includes inactive lanes. Integer-only bounds reject zero, subnormal, Inf/NaN and excessive operands before any new FP operation; every other vector retains original execution. No active masks, operand loading, aliasing writes, identity tests, metadata or caller operations change.'),
        reference=dict(url='https://www.tuhh.de/ti3/paper/rump/OgRuOi05.pdf',
                       location='Algorithms3.2/3.3 and Theorem3.4',
                       role='Original split/TwoProduct exact identity when product underflow is excluded; bounds and binary64 equivalence above are this source derivation.'),
        performance_scope='The fast branch replaces17 binary32 arithmetic operations with one original binary32 multiply, two binary64 operations and four exact vector conversions, plus integer bounds/branches. It reaches three product32 calls per general multiply and the original repeated multiply calls inside divide/sqrt. Actual whole-eight/nonzero-residual hit rate, code/register/stack effects and benefit are unmeasured; mixed exceptional/zero vectors retain the old chain. No claim that this alone closes the~3.5x frame gap.',
        scope='One internal primitive only in a separately derived optional-v4 source. Original retained8.cpp, original-v3 object, baseline dispatcher, five-word ABI, eight32-bit lanes, all operations outside product32, refinements, masks, error allowance, raw3/scalar/provider/host/SPIRV and acceptance criteria are unchanged. No compiler execution, link, test, profiling, provider mutation or qualification claim.')
    with (W / 'derivation.json').open('x') as stream:
        json.dump(proof, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(source_only=True, inverse_exact=True, changed_function='product32',
                          compiled=False, executed=False)))


if __name__ == '__main__':
    assert sys.argv[1:] == ['derive']
    derive()
