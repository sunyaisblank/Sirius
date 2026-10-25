"""Lower three pure retained primitives to lane selections without new arithmetic."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import sys

R = Path.cwd(); W = Path(__file__).resolve().parent
assert sys.argv[1:] in ([], ['--write'])
pin = lambda data: {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write']
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)

def function(text, name):
    match = re.search(r'^static (?:inline )?[^\n]+\b' + name + r'\([^\n]+\)\n\{', text, re.M)
    assert match, name
    depth = 1; end = match.end()
    while depth:
        depth += (text[end] == '{') - (text[end] == '}'); end += 1
    return text[match.start():end]

def lower(source):
    changes = {}
    name = 'RetainedSum32_0'; old = function(source, name)
    args = re.search(r'\(float (\w+), float (\w+), varying float \* uniform (\w+), varying float \* uniform (\w+)\)', old).groups()
    a, b, hi, lo = args
    values = re.search(r'float (_S\d+);[\s\S]*?float (_S\d+);', old).groups()
    left, right = values
    ordered = re.search(r'    float (_S\d+) = OrderedAdd32_0\(' + left + ', ' + right + r'\);', old)
    assert ordered
    # Keep the exact original three ordered floating expressions and stores.
    arithmetic = old[ordered.start():old.rindex('    return;')]
    new = old[:old.index('{') + 1] + f'''
    bool SiriusSwap = ((RPBits_0({a})) & 2147483647U) < ((RPBits_0({b})) & 2147483647U);
    float {left} = select(SiriusSwap, {b}, {a});
    float {right} = select(SiriusSwap, {a}, {b});
''' + arithmetic + '    return;\n}'
    assert arithmetic in old and arithmetic in new
    changes[name] = (old, new)

    name = 'RPUpAdd_0'; old = function(source, name)
    a, b = re.search(r'\(float (\w+), float (\w+)\)', old).groups()
    expression = re.search(r'return (RPFromBits_0\(RPBits_0\(OrderedAdd32_0\([^;]+);', old)[1]
    new = old[:old.index('{') + 1] + f'''
    bool SiriusBounded = RPBounded_0({a}, 2105540608U) & RPBounded_0({b}, 2105540608U);
    float SiriusRounded = {expression};
    float SiriusNonzeroA = select(RPScalarEqual_0({b}, RPZero_0()), {a}, SiriusRounded);
    float SiriusSelected = select(RPScalarEqual_0({a}, RPZero_0()), {b}, SiriusNonzeroA);
    return select(SiriusBounded, SiriusSelected, RPFailureBound_0());
}}'''
    assert expression in old and expression in new
    changes[name] = (old, new)

    name = 'RTAccumulateTerm_0'; old = function(source, name)
    a, x = re.search(r'\(varying RetainedTriple_0 \* uniform (\w+), float (\w+)\)', old).groups()
    sums = re.findall(r'RetainedSum32_0\([^;]+, &(\w+), &(\w+)\);', old)
    assert len(sums) == 3
    hi, lo, tail = [v[0] for v in sums]
    error = re.search(r'    float (\w+) = RPUpAdd_0\([^;]+;', old)
    assert error
    error_name = error[1]
    arithmetic = old[old.index('    float ' + hi + ';'):error.end()]
    assert len(re.findall(r'RetainedSum32_0\(', arithmetic)) == 3
    zero_fields = '\n'.join(f'''    SiriusZero.{field} = select(SiriusPositiveZero & ((RPBits_0({a}->{field}) & 2147483647U) == 0U), RPZero_0(), {a}->{field});''' for field in ('hi_0', 'lo_0', 'tail_0'))
    input_bounds = ' &\n        '.join(f'RPBounded_0({a}->{field}, 2105540608U)' for field in ('hi_0', 'lo_0', 'tail_0', 'error_0'))
    output_bounds = ' &\n        '.join(f'RPBounded_0({value}, 2105540608U)' for value in (hi, lo, tail, error_name))
    result_fields = '\n'.join(f'''    SiriusResult.{field} = select(SiriusGood, select(SiriusIsZero, SiriusZero.{field}, {value}), RPZero_0());''' for field, value in zip(('hi_0', 'lo_0', 'tail_0', 'error_0'), (hi, lo, tail, error_name)))
    new = old[:old.index('{') + 1] + f'''
    bool SiriusIsZero = (RPBits_0({x}) & 2147483647U) == 0U;
    bool SiriusZeroGood = ({a}->valid_0 != 0U) &
        {input_bounds};
    bool SiriusPositiveZero = RPBits_0({x}) == 0U;
    RetainedTriple_0 SiriusZero;
{zero_fields}
    SiriusZero.error_0 = select((RPBits_0({a}->error_0) & 2147483647U) == 0U, RPZero_0(), {a}->error_0);
    SiriusZero.valid_0 = 1U;
    // Preserve the cheap zero result when every currently executing lane has it.
    if (all(SiriusIsZero & SiriusZeroGood))
        return SiriusZero;

''' + arithmetic + f'''
    bool SiriusGeneralGood = ({a}->valid_0 != 0U) & RPBounded_0({x}, 2071986176U) &
        {output_bounds};
    bool SiriusGood = (SiriusIsZero & SiriusZeroGood) | ((!SiriusIsZero) & SiriusGeneralGood);
    RetainedTriple_0 SiriusResult;
{result_fields}
    SiriusResult.valid_0 = select(SiriusGood, 1U, 0U);
    return SiriusResult;
}}'''
    assert arithmetic in old and arithmetic in new
    changes[name] = (old, new)
    candidate = source
    for old, new in changes.values():
        assert candidate.count(old) == 1; candidate = candidate.replace(old, new, 1)
    inverse = candidate
    for old, new in changes.values():
        assert inverse.count(new) == 1; inverse = inverse.replace(new, old, 1)
    assert inverse == source
    return candidate, changes

basis_path = R / 'out/native-cpu-product-inline/ispc-derivation.json'
basis = json.loads(basis_path.read_text())
inputs = {str(basis_path.relative_to(R)): pin(basis_path.read_bytes()),
          str(Path(__file__).resolve().relative_to(R)): pin(Path(__file__).read_bytes())}
mesa = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
for name in ('src/util/u_math.h', 'src/gallium/drivers/llvmpipe/lp_state_cs.c'):
    p = mesa / name; inputs[str(p.relative_to(R))] = pin(p.read_bytes())
assert 'util_fpstate_set((saved_fpstate | 0x1f80u) & ~0xe040u);' in (mesa / 'src/util/u_math.h').read_text()
worker = function((mesa / 'src/gallium/drivers/llvmpipe/lp_state_cs.c').read_text().replace('static void\ncs_exec_fn(', 'static void cs_exec_fn('), 'cs_exec_fn')
assert worker.index('util_fpstate_get()') < worker.index('util_fpstate_set_ieee32(fpstate)') < worker.index('lp_native_cpu_packet(') < worker.index('util_fpstate_set(fpstate)')
stages = {}; patch = ''; count = 0
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']; data = path.read_bytes(); assert pin(data) == entry['derived']
    source = data.decode(); candidate, changes = lower(source)
    literal = lambda s: [line for line in s.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert literal(source) == literal(candidate)
    count += len(literal(source))
    destination = R / 'bin/linux-gcc/native-cpu-predicated-arithmetic/ispc-source' / (slug + '.ispc')
    publish(destination, candidate.encode())
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {**entry, 'predicated_basis_path': str(path.relative_to(R)), 'predicated_basis': pin(data),
        'derived_path': str(destination.relative_to(R)), 'derived': pin(candidate.encode()),
        'modified_functions': {name: {'basis': pin(old.encode()), 'candidate': pin(new.encode())} for name, (old, new) in changes.items()},
        'inverse_restores_exact_W_source': True, 'literal_instruction_calls_preserved': len(literal(source))}
    patch += ''.join(difflib.unified_diff(source.splitlines(True), candidate.splitlines(True), fromfile='a/' + slug + '.ispc', tofile='b/' + slug + '.ispc'))
assert count == 9380
result = {'instructions': count, 'stages': stages, 'source_inputs': inputs,
    'scope': 'Only RetainedSum32 operand sorting, RPUpAdd result choice and RTAccumulateTerm lane result/validity selection change. Original ordered floating expressions and their live-result order retained; no arithmetic rederivation, FP/ISA/eligibility/count/field/column/RHS/trial/quartic/output/admission changes. U scalar/F/I/Q/V deferred-NIR paths remain. Source reasoning, not machine-checked equivalence or runtime qualification.',
    'reasoning': [
        'RetainedSum chooses the identical two operands by the same unsigned magnitude-bit predicate; ordered add/subtract expressions and stores are copied byte-for-byte.',
        'RPUpAdd keeps bounds rejection first in the selected result, then a==zero before b==zero (preserving signed zero), then the exact original ordered-add bit increment. Bitwise select preserves chosen payload bits.',
        'Term retains the exact S signed-zero result and bounds, including canonical positive error zero and valid=1. The uniform all-active-lanes-zero branch returns it cheaply. Other lanes run the original three ordered sums and upward error addition in the original order; bitwise validity selects exactly the old generic result or all-zero invalid tuple.',
        'Eager expressions in unselected/refused lanes contain only local bit/arithmetic operations, no indexing, division, libm/errno or external state access. They may set extra masked SSE status flags. Existing cs_exec_fn masks all exceptions, establishes RTE/FTZ/DAZ-off and restores the full saved FP state for both normal and refused arithmetic. No shader primitive observes sticky flags; status changes cannot escape the owned scope.',
        'No eligibility precondition is strengthened. Arbitrary divergent lane values retain the generic selected result/refusal, and inactive masks remain governed by original export/caller writes. Full38/originalC8 and later acceptance remain required.'
    ],
    'ISPC_selection_reference': 'https://ispc.github.io/ispc.html#logical-and-selection-operations'}
publish(W / 'arithmetic.patch', patch.encode())
publish(W / 'ispc-derivation.json', (json.dumps(result, indent=2) + '\n').encode())
print(json.dumps({'complete_stages': 3, 'modified_primitives_each': 3, 'literal_instructions': count, 'runtime_qualified': False}))
