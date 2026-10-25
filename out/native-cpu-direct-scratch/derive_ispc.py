"""Prepare a contingent result-to-scratch ABI; never compile or run a workload."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import sys

R = Path.cwd(); W = Path(__file__).resolve().parent
assert sys.argv[1:] in ([], ['--write'])
pin = lambda b: {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}

def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write']
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)

def function(text, name):
    m = re.search(r'^static (?:inline )?[^\n]+\b' + name + r'\([^\n]+\)\n\{', text, re.M); assert m, name
    depth = 1; end = m.end()
    while depth:
        depth += (text[end] == '{') - (text[end] == '}'); end += 1
    return text[m.start():end]

basis_path = R / 'out/native-cpu-predicated-arithmetic/ispc-derivation.json'
basis = json.loads(basis_path.read_text())
inputs = {str(basis_path.relative_to(R)): pin(basis_path.read_bytes()),
          str(Path(__file__).resolve().relative_to(R)): pin(Path(__file__).read_bytes())}
stages = {}; patch = ''; count = 0
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']; data = path.read_bytes(); assert pin(data) == entry['derived']
    source = data.decode(); store_helpers = []; helper_proofs = {}
    read = function(source, 'ReadValueUniform_0')
    read_parameters = re.search(r'ReadValueUniform_0\(uniform uint32_t (\w+), uniform uint32_t (\w+), varying KernelContext_0 \* uniform (\w+)\)', read).groups()
    for name in ('RTAdd_0', 'RTMultiply_0'):
        original = function(source, name)
        a, b = re.search(r'\(varying RetainedTriple_0 \* uniform (\w+), varying RetainedTriple_0 \* uniform (\w+)\)', original).groups()
        body = original[original.index('{') + 1:-1]
        returns = list(re.finditer(r'\breturn\s+([^;]+);', body)); assert len(returns) >= 2
        transformed = body
        inverse_replacements = []
        for i, match in reversed(list(enumerate(returns))):
            replacement = (f'{{ RetainedTriple_0 SiriusStoredResult_{i} = ' + match[1] + ';\n'
                f'        WriteValueUniform_0(SiriusDestination, &SiriusStoredResult_{i}, SiriusRowBase, SiriusContext);\n'
                '        return; }')
            transformed = transformed[:match.start()] + replacement + transformed[match.end():]
            inverse_replacements.append((replacement, match[0]))
        inverse = transformed
        for replacement, old in inverse_replacements:
            assert inverse.count(replacement) == 1; inverse = inverse.replace(replacement, old, 1)
        assert inverse == body
        new_name = name.replace('_0', 'Store_0')
        helper = f'''static noinline void {new_name}(uniform uint32_t SiriusDestination,
    uniform uint32_t SiriusA, uniform uint32_t SiriusB, uniform uint32_t SiriusRowBase,
    varying KernelContext_0 * uniform SiriusContext)
{{
    // Load both full operands before any destination write, including aliases.
    RetainedTriple_0 SiriusReadA = ReadValueUniform_0(SiriusA, SiriusRowBase, SiriusContext);
    RetainedTriple_0 SiriusReadB = ReadValueUniform_0(SiriusB, SiriusRowBase, SiriusContext);
    varying RetainedTriple_0 * uniform {a} = &SiriusReadA;
    varying RetainedTriple_0 * uniform {b} = &SiriusReadB;
''' + transformed + '\n}\n'
        store_helpers.append(helper)
        helper_proofs[name] = {'basis': pin(original.encode()), 'store_helper': pin(helper.encode()),
            'original_return_expressions': [m[1] for m in returns], 'inverse_restores_exact_arithmetic_body': True,
            'all_result_words_use_original_WriteValueUniform': True}
    added = '\n// Direct private-scratch result ABI; original arithmetic helpers retained.\n' + '\n'.join(store_helpers) + '\n'
    write = function(source, 'WriteValueUniform_0'); position = source.index(write) + len(write)
    candidate = source[:position] + added + source[position:]
    evaluate = function(candidate, 'EvaluateOneStatic_0')
    d, a, b, c, row, rowbase, regs, ctx = re.search(r'EvaluateOneStatic_0\(uniform uint32_t (\w+), uniform uint32_t (\w+), uniform uint32_t (\w+), uniform uint32_t (\w+), uint32_t (\w+), uniform uint32_t (\w+), uint32_t (\w+), varying KernelContext_0 \* uniform (\w+)\)', evaluate).groups()
    marker = re.search(r'^    RetainedTriple_0 \w+ = RTInvalid_0\(\);', evaluate, re.M); assert marker
    fast = f'''    if (StaticOp == 2U || StaticOp == 4U) {{
        if ({a} >= {regs} || {b} >= {regs})
            return false;
        if (StaticOp == 2U)
            RTAddStore_0({d}, {a}, {b}, {rowbase}, {ctx});
        else
            RTMultiplyStore_0({d}, {a}, {b}, {rowbase}, {ctx});
        return true;
    }}

'''
    altered = evaluate[:marker.start()] + fast + evaluate[marker.start():]
    assert candidate.count(evaluate) == 1; candidate = candidate.replace(evaluate, altered, 1)
    inverse = candidate.replace(altered, evaluate, 1).replace(added, '', 1)
    assert inverse == source
    literal = lambda s: [line for line in s.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert literal(candidate) == literal(source)
    count += len(literal(source))
    destination = R / 'bin/linux-gcc/native-cpu-direct-scratch/ispc-source' / (slug + '.ispc')
    publish(destination, candidate.encode())
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {**entry, 'direct_scratch_basis_path': str(path.relative_to(R)), 'direct_scratch_basis': pin(data),
        'derived_path': str(destination.relative_to(R)), 'derived': pin(candidate.encode()), 'helpers': helper_proofs,
        'literal_instruction_calls_preserved': len(literal(source)), 'inverse_restores_exact_X_source': True}
    patch += ''.join(difflib.unified_diff(source.splitlines(True), candidate.splitlines(True), fromfile='a/' + slug + '.ispc', tofile='b/' + slug + '.ispc'))
assert count == 9380
proof = {'instructions': count, 'stages': stages, 'source_inputs': inputs,
    'status': 'contingent source preparation only; no compilation/workload/qualification',
    'scope': 'Preserve all literal tuples and original Add/Multiply arithmetic bodies/branch order/return expressions. New shared noinline void helpers load both full operands before any write and materialize each returned tuple through original WriteValueUniform directly into its private scratch destination. This avoids a caller-side returned-struct temporary without broad whole-operator inlining. Original helpers remain for other consumers. No count/field/precision/admission change.',
    'reasoning': [
        'Original destination/opcode guards remain first. New specialization for StaticOp2/4 retains a/b register bounds before writes. Canonical full admission proves the static metadata domain; other opcodes keep their original body.',
        'All five words of both operands are read before destination writes, so reused/aliased canonical registers observe the same prior values. Every original return expression is copied into a local and sent through the original five-field scratch writer; an exact inverse restores the arithmetic body.',
        'Earlier writes for already-returned lanes touch only their private scratch components. Other lanes use operands captured before those writes; no cross-lane scratch dependency or external publication is introduced. Original export publication/refusal/whole-state rollback remains. Compilation/code and original numerical/runtime gates are still required.',
        'Do not compile or execute this contingent candidate unless Y fails its decisive original frame. Remove unused generated sources if Y succeeds.'
    ]}
publish(W / 'direct-scratch.patch', patch.encode())
publish(W / 'ispc-derivation.json', (json.dumps(proof, indent=2) + '\n').encode())
print(json.dumps({'complete_source_copies': 3, 'shared_store_helpers_each': 2, 'literal_instructions': count, 'compiled': False}))
