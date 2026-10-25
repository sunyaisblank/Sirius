"""Prepare an unapplied exact signed-zero accumulation specialization for all3 stages."""
from pathlib import Path
import difflib
import hashlib
import json
import re

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'out/native-cpu-accumulate-inline'


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(name, value):
    path = W / name
    assert not path.exists()
    path.write_text(value)


derived = json.loads((B / 'ispc-derivation.json').read_text())
compiled = json.loads((B / 'ispc-inline-compile-result.json').read_text())
assert derived['instructions'] == 9380 and compiled['all_objects_compiled']
assert json.loads((B / 'inline-code-inspection.json').read_text())['passed']
patch, stages, calls = '', {}, 0
inputs = {}
for name in ('ispc-derivation.json', 'ispc-inline-compile-result.json', 'inline-code-inspection.json'):
    path = B / name
    inputs[str(path.relative_to(R))] = pin(path.read_bytes())
for slug in ('transport', 'endpoint', 'dopri_phase'):
    entry = derived['stages'][slug]
    path = R / entry['derived_path']
    data = path.read_bytes()
    assert pin(data) == entry['derived']
    original = data.decode()
    header = re.search(r'static inline RetainedTriple_0 RTAccumulateTerm_0\(varying RetainedTriple_0 \* uniform (\w+), float (\w+)\)', original)
    assert header
    a, x = header.groups()
    first_return = original.index('        return RTInvalid_0();\n    }', header.end())
    position = first_return + len('        return RTInvalid_0();\n    }')
    assert 'RPBounded_0(' + x + ', 2071986176U)' in original[header.end():position]
    assert '->valid_0' in original[header.end():position]
    assert 'SiriusZeroTerm' not in original
    arm = f'''

    // Exact original outputs for either sign of zero, under the existing IEEE32 scope.
    if ((RPBits_0({x}) & 2147483647U) == 0U)
    {{
        if (!RPBounded_0({a}->hi_0, 2105540608U) ||
            !RPBounded_0({a}->lo_0, 2105540608U) ||
            !RPBounded_0({a}->tail_0, 2105540608U) ||
            !RPBounded_0({a}->error_0, 2105540608U))
            return RTInvalid_0();
        bool SiriusZeroTermPositive = RPBits_0({x}) == 0U;
        RetainedTriple_0 SiriusZeroTermResult;
'''
    for field in ('hi_0', 'lo_0', 'tail_0'):
        arm += f'        SiriusZeroTermResult.{field} = (SiriusZeroTermPositive && ((RPBits_0({a}->{field}) & 2147483647U) == 0U)) ? RPZero_0() : {a}->{field};\n'
    arm += f'''        SiriusZeroTermResult.error_0 = ((RPBits_0({a}->error_0) & 2147483647U) == 0U) ? RPZero_0() : {a}->error_0;
        SiriusZeroTermResult.valid_0 = 1U;
        return SiriusZeroTermResult;
    }}'''
    candidate = original[:position] + arm + original[position:]
    assert candidate.replace(arm, '', 1) == original and candidate.count(arm) == 1
    literal = lambda text: [line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert literal(candidate) == literal(original)
    assert len(literal(original)) == entry['literal_instruction_calls_preserved']
    calls += len(literal(original))
    future = 'bin/linux-gcc/native-cpu-zero-term/ispc-source/' + slug + '.ispc'
    patch += ''.join(difflib.unified_diff(original.splitlines(True), candidate.splitlines(True),
        fromfile='a/' + future, tofile='b/' + future, n=3))
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {'basis_path': str(path.relative_to(R)), 'basis': pin(data),
                    'future_candidate_path': future, 'candidate': pin(candidate.encode()),
                    'inverse_removes_one_zero_arm_only': True,
                    'literal_instruction_calls_preserved': len(literal(original)),
                    'zero_arm': arm}
assert calls == 9380
proof = {'status': 'unapplied source proposal; no new source copies, compilation, link, workload or qualification',
         'instructions': calls, 'stages': stages, 'source_inputs': inputs,
         'premise': 'R removed Term ABI boundaries but retained all three sorted sums and radius/bounds work for a zero term. RTAdd of a nonzero RTConstant still passes its exactly-zero low/tail components through two such accumulations. Source-supported opportunity; dynamic frequency, compiler code shape and benefit remain unmeasured. No parameter/flag sweep.',
         'scope': 'All3 full stage programs retain every literal instruction call and all fields/columns/RHS/RK/quartic/guard/layout/counts/FP/eligibility. Original prefix rejects invalid a.valid and out-of-bound x exactly as before. For nonzero x, every original body byte remains and executes. No unmasked qualifier, precision/budget/cache/count/reference/governor change.',
         'semantic_argument': [
             'Existing scope is binary32 RTE, FTZ/DAZ off, no contraction. For x=+0 or -0 and any finite a, sorted RetainedSum32 selects a first because absbits(a)<absbits(x) is false. s=a+x; s-a is exactly +0 in every sign-zero case and for every nonzero finite a, including subnormals. Its residual x-(s-a) is therefore exactly x, with its sign.',
             'Induct through all3 sums: the same signed zero is carried through hi/lo/tail. For x=+0 each center component keeps its bits except -0 becomes +0. For x=-0 all3 center component bits remain unchanged. No component can cross an original magnitude bound by adding zero.',
             'The final residual absolute value is +0. RPUpAdd(error,+0) first applies the unchanged magnitude bound; if error is either zero it returns the second operand +0, otherwise returns error unchanged. Thus the error field only canonicalizes -0 to +0.',
             'If any hi/lo/tail/error magnitude exceeds the original 2105540608 bound, the original result is canonical RTInvalid; the zero arm returns exactly that same value. This also covers nonfinite and signaling NaN inputs as output-bit semantics. If a.valid is nonzero, the original successful result always writes valid=1, and the zero arm does too.',
             'No external FP-environment equivalence is claimed outside the existing native scoped IEEE32 call, which restores the complete caller MXCSR. Internal sticky exception flags may differ for rejected nonfinite inputs; they are not exposed by that scope. Exact output bits, radius/enclosures, refusal and active/inactive ownership remain required.'
         ],
         'signed_zero_table': [
             {'term': '+0', 'input_component': '+0', 'sum': '+0', 'residual': '+0'},
             {'term': '+0', 'input_component': '-0', 'sum': '+0', 'residual': '+0'},
             {'term': '-0', 'input_component': '+0', 'sum': '+0', 'residual': '-0'},
             {'term': '-0', 'input_component': '-0', 'sum': '-0', 'residual': '-0'}],
         'required_before_adoption': 'Review the all-binary32 bit/rounding argument; actual all3 same-tool/flags compilation and ISA/code checks; fresh changed-Transport all38/all187960-word compatibility, original C8/all32/51markers/54joins, sole ABBA22+Decimal, all8 controls and original300s frame. Reuse the exact existing driver/baseline diagnostic objects if their source/ABI/input/tool provenance still matches; no reason to compile their unchanged C++ again. Full science/release/main qualification remains separate.'}
publish('zero-term.patch', patch)
publish('source-proof.json', json.dumps(proof, indent=2) + '\n')
print(json.dumps({'source_only': True, 'stages': 3, 'instructions': calls, 'patch_applied': False, 'workloads_executed': 0}))
