"""Emit integer-only value16 adapters to the exact existing strict value8 cores.

This writes source only. It never compiles, loads a module, or runs arithmetic.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re

R = Path(__file__).resolve().parents[2]
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-wide-private'
MATH = R / 'bin/linux-gcc/native-cpu-local-effects'
PINS = {
    'math.bc': (296980, '93e076866edeb28c329a65a81dad1485a2717a8104ec4b7879002c78a9e6a6a1'),
    'math.ll': (2789102, '8c41b14d6e30a9befc48be12c59324840421f51557c245595575ec4105631ed0'),
    'math.o': (259216, '69f9d771c5cf2d0de63e32ac70ac534961e9283a02753047784902e3ee54cdbe'),
}
RECIPES = ('dense_add', 'dense_multiply', 'dense_divide', 'dense_sqrt',
           'tagged_add', 'tagged_subtract', 'tagged_multiply', 'tagged_divide', 'tagged_sqrt')
V8, V16 = '<8 x i32>', '<16 x i32>'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def vector(indices):
    return '<' + ', '.join('i32 ' + str(i) for i in indices) + '>'


def aggregate(fields, lane_type):
    return '{ ' + ', '.join([lane_type] * fields) + ' }'


def recipe(name):
    tagged = name.startswith('tagged_')
    fields, arguments = (6, 12) if tagged else (5, 10)
    small, large = aggregate(fields, V8), aggregate(fields, V16)
    old = 'sirius_local_value_' + name + '8'
    new = 'sirius_wide_value_' + name + '16'
    old_params = ', '.join([V8] * arguments + ['i32'])
    new_params = ', '.join([f'{V16} %a{i}' for i in range(arguments)] + ['i32 %mask'])
    lines = [f'declare hidden x86_regcallcc {small} @{old}({old_params}) #1', '',
             f'define hidden x86_regcallcc {large} @{new}({new_params}) #0 {{',
             'entry:', '  %mask.lo = and i32 %mask, 255',
             '  %mask.shift = lshr i32 %mask, 8', '  %mask.hi = and i32 %mask.shift, 255',
             '  %has.lo = icmp ne i32 %mask.lo, 0',
             '  br i1 %has.lo, label %call.lo, label %join.lo', 'call.lo:']
    for half, indices in (('lo', range(8)), ('hi', range(8, 16))):
        if half == 'hi':
            lines += ['call.hi:']
        for i in range(arguments):
            lines += [f'  %{half}.a{i} = shufflevector {V16} %a{i}, {V16} zeroinitializer, {V8} {vector(indices)}']
        args = ', '.join([f'{V8} %{half}.a{i}' for i in range(arguments)] + [f'i32 %mask.{half}'])
        lines += [f'  %{half}.result = call x86_regcallcc {small} @{old}({args}) #1',
                  f'  br label %join.{half}', f'join.{half}:']
        previous = 'entry' if half == 'lo' else 'join.lo'
        lines += [f'  %{half}.value = phi {small} [ zeroinitializer, %{previous} ], [ %{half}.result, %call.{half} ]']
        if half == 'lo':
            lines += ['  %has.hi = icmp ne i32 %mask.hi, 0',
                      '  br i1 %has.hi, label %call.hi, label %join.hi']
    for i in range(fields):
        lines += [f'  %lo.field{i} = extractvalue {small} %lo.value, {i}',
                  f'  %hi.field{i} = extractvalue {small} %hi.value, {i}',
                  f'  %field{i} = shufflevector {V8} %lo.field{i}, {V8} %hi.field{i}, {V16} {vector(range(16))}']
        prior = 'zeroinitializer' if i == 0 else f'%result{i-1}'
        lines += [f'  %result{i} = insertvalue {large} {prior}, {V16} %field{i}, {i}']
    lines += [f'  ret {large} %result{fields-1}', '}', '']
    body = '\n'.join(lines)
    assert body.count(' = call x86_regcallcc ') == 2
    assert not re.search(r'\b(?:undef|poison|freeze|fadd|fsub|fmul|fdiv|alloca|load|store)\b', body)
    return body, dict(old=old, new=new, vector_arguments=arguments, return_vectors=fields,
                      low_mask='mask & 255', high_mask='(mask >> 8) & 255',
                      zero_half='Skip the call and retain the original all-zero private destination.',
                      order='low call, then high call; each only when its mask is nonzero')


def sources():
    for name, (size, digest) in PINS.items():
        assert pin(MATH / name) == dict(bytes=size, sha256=digest), name
    text = (MATH / 'math.ll').read_text()
    layout, = re.findall(r'^target datalayout = (.+)$', text, re.M)
    triple, = re.findall(r'^target triple = (.+)$', text, re.M)
    expected = []
    for name in RECIPES:
        fields, arguments = (6, 12) if name.startswith('tagged_') else (5, 10)
        prefix = f'define hidden x86_regcallcc {aggregate(fields, V8)} @sirius_local_value_{name}8('
        header, = [line for line in text.splitlines() if line.startswith(prefix)]
        params = header[len(prefix):header.index(') local_unnamed_addr')]
        assert params == ', '.join([f'{V8} %{i}' for i in range(arguments)] + [f'i32 %{arguments}'])
        expected.append(header)
    assert len(re.findall(r'^define .*@sirius_local_value_', text, re.M)) == 9
    bodies, proofs = zip(*(recipe(name) for name in RECIPES))
    source = ('; Private SIMD16 interface; all arithmetic remains in the original strict8 object.\n'
              f'target datalayout = {layout}\ntarget triple = {triple}\n\n' + '\n'.join(bodies) +
              '\nattributes #0 = { noinline nounwind strictfp memory(inaccessiblemem: readwrite) uwtable '
              '"target-cpu"="x86-64-v4" "min-legal-vector-width"="512" "prefer-vector-width"="512" }\n'
              'attributes #1 = { noinline nounwind strictfp memory(inaccessiblemem: readwrite) }\n')
    assert source.count('define hidden ') == source.count('declare hidden ') == 9
    assert source.count(' = call x86_regcallcc ') == 18
    contract = dict(input_math={name:pin(MATH / name) for name in ('math.bc', 'math.ll', 'math.o')},
        exact_old_core_headers=expected, recipes=list(proofs),
        initialized_return='Both phi zero paths and every joined result field are fully defined; no undef, poison or freeze.',
        scope='Private guarded positive-job SIMD16 only, within original complete saved/masked/restored FP scope. '
              'Original strict eight-lane math object, arithmetic/rounding/fast and slow decisions, public math and fallbacks stay exact. '
              'This adapter has no FP instructions or ordinary memory effects. Raw16 source/mask/subgroup/context/native/admission proof remains required.',
        qualification='Source only; no LLVM parse, compiler, native ABI, numerical, timing, provider or frame qualification.')
    return source, contract


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('check', 'derive'))
    args = parser.parse_args()
    source, contract = sources()
    if args.phase == 'check':
        print(json.dumps(dict(source_only=True, wrappers=9, guarded_original_calls=18,
                              source_bytes=len(source.encode()), source_sha256=hashlib.sha256(source.encode()).hexdigest())))
        return
    output = B / 'wide-math.ll'
    receipt = W / 'math-source-contract.json'
    assert not output.exists() and not receipt.exists()
    B.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        stream.write(source)
    contract['output'] = dict(path=str(output.relative_to(R)), **pin(output))
    with receipt.open('x') as stream:
        json.dump(contract, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(source_only=True, output=contract['output'])))


if __name__ == '__main__':
    main()
