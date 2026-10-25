"""Prepare constant canonical register metadata; no compilation or workload."""
from pathlib import Path
import hashlib
import json
import re
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
assert sys.argv[1:] in ([], ['--write'])


def identity(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write']
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)


def function(source, name):
    match = re.search(r'^static [^\n]*\b' + name + r'\(', source, re.M)
    assert match, name
    opening = source.index('{', match.end())
    end, depth = opening + 1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[match.start():end]


basis_path = R / 'out/native-cpu-direct-scratch/ispc-derivation.json'
basis = json.loads(basis_path.read_text())
inputs = {str(basis_path.relative_to(R)): identity(basis_path.read_bytes()),
          str(Path(__file__).resolve().relative_to(R)): identity(Path(__file__).read_bytes())}
stages, instruction_count, layer_count = {}, 0, 0
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']
    data = path.read_bytes()
    assert identity(data) == entry['derived']
    original = data.decode()
    candidate = original
    changes = []

    def replace(before, after):
        global candidate
        assert candidate.count(before) == 1
        candidate = candidate.replace(before, after, 1)
        changes.append((before, after))

    known = function(original, 'EvaluateKnownLayer_0')
    headers = list(re.finditer(
        r'if\(programBase == tableBase \+ (\d+)U && count == (\d+)U && registers == (\d+)U\)', known))
    assert len(headers) in (1, 2)
    families = {}
    for i, header in enumerate(headers):
        branch = known[header.end():headers[i + 1].start() if i + 1 < len(headers) else len(known)]
        calls = re.findall(r'EvaluateStaticLayer_(\d+)_(\d+)_0\(row, rowBase, registers, context\);', branch)
        assert calls
        family, = {int(call[0]) for call in calls}
        assert family not in families
        families[family] = {'program_offset': int(header[1]), 'instructions': int(header[2]),
                            'registers': int(header[3]), 'layers': [int(call[1]) for call in calls]}

    replace('template<uint32_t StaticOp>\nstatic bool EvaluateOneStatic_0(',
            'template<uint32_t StaticOp, uint32_t StaticRegisters>\nstatic bool EvaluateOneStatic_0(')
    one = function(candidate, 'EvaluateOneStatic_0')
    parameter, = re.findall(r', uint32_t (\w+), varying KernelContext_0 \* uniform \w+\)', one.split('{', 1)[0])
    specialized = one.replace(', uint32_t ' + parameter + ', varying KernelContext_0',
                              ', varying KernelContext_0', 1)
    specialized = specialized.replace('{', '{\n    uniform const uint32_t ' + parameter + ' = StaticRegisters;', 1)
    replace(one, specialized)

    tuples_before, tuples_after = [], []
    for family, metadata in families.items():
        for layer in metadata['layers']:
            name = f'EvaluateStaticLayer_{family}_{layer}_0'
            assert len(re.findall(r'\b' + name + r'\(', original)) == 2
            assert len(re.findall(r'\b' + name + r'\(', known)) == 1
            before = function(candidate, name)
            signature = before.split('{', 1)[0]
            registers, = re.findall(r', uint32_t (\w+), varying KernelContext_0 \* uniform \w+\)', signature)
            after = before.replace(', uint32_t ' + registers + ', varying KernelContext_0',
                                   ', varying KernelContext_0', 1)
            pattern = (r'EvaluateOneStatic_0<(\d+)U>\((\d+)U, (\d+)U, (\d+)U, (\d+)U, '
                       r'(\w+), (\w+), ' + re.escape(registers) + r', (\w+)\)')

            def instruction(match):
                op, dst, a, b, c = map(int, match.groups()[:5])
                assert 0 <= op <= 11 and dst < metadata['registers']
                arity = {0: 0, 1: 0, 2: 2, 3: 2, 4: 2, 5: 2,
                         6: 1, 7: 1, 8: 1, 9: 3, 10: 3, 11: 3}[op]
                assert all(index < metadata['registers'] for index in (a, b, c)[:arity])
                tuples_before.append((family, layer, op, dst, a, b, c))
                row, base, context = match.groups()[5:]
                result = (f'EvaluateOneStatic_0<{op}U, {metadata["registers"]}U>'
                          f'({dst}U, {a}U, {b}U, {c}U, {row}, {base}, {context})')
                tuples_after.append((family, layer, op, dst, a, b, c))
                return result

            after, number = re.subn(pattern, instruction, after)
            assert number > 0 and 'registers' not in after
            replace(before, after)
            layer_count += 1

    changed_known = re.sub(r'(EvaluateStaticLayer_\d+_\d+_0\(row, rowBase), registers(, context\);)',
                           r'\1\2', known)
    assert changed_known != known
    replace(known, changed_known)
    inverse = candidate
    for before, after in reversed(changes):
        assert inverse.count(after) == 1
        inverse = inverse.replace(after, before, 1)
    assert inverse == original
    assert tuples_before == tuples_after and len(tuples_before) == sum(v['instructions'] for v in families.values())
    # Structural fallback, numerical helpers, row lifecycle and publication remain byte-identical.
    for name in ('EvaluateOne_0', 'EvaluateProgram_0', 'ComputeOwnedRow_0', 'RTAddStore_0', 'RTMultiplyStore_0'):
        assert function(candidate, name) == function(original, name)
    destination = R / 'bin/linux-gcc/native-cpu-static-registers/ispc-source' / (slug + '.ispc')
    publish(destination, candidate.encode())
    inputs[str(path.relative_to(R))] = identity(data)
    stages[slug] = {'basis_path': str(path.relative_to(R)), 'basis': identity(data),
                    'derived_path': str(destination.relative_to(R)), 'derived': identity(candidate.encode()),
                    'families': families, 'instructions': len(tuples_before),
                    'inverse_restores_exact_Z_source': True, 'ordered_instruction_tuples_unchanged': True}
    instruction_count += len(tuples_before)
assert instruction_count == 9380 and layer_count == 399
proof = {'status': 'contingent source only; not compiled or qualified', 'instructions': instruction_count,
         'layers': layer_count, 'stages': stages, 'source_inputs': inputs,
         'scope': 'Canonical per-family register counts become compile-time uniform constants in the shared opcode bodies. Original exact family/header/layer guards still dominate every static layer call. Generic fallback, all arithmetic bodies, ordered literal tuples, numerical refusal, outputs, fields, trials and publication are unchanged; exact whole-source inverse checked. No runtime input eligibility is strengthened.',
         'reasoning': [
             'The original EvaluateKnownLayer branch requires the actual register count to equal the canonical family constant before any layer call. Layer helpers are static and referenced only by those admitted branches. Replacing the propagated varying count with that exact constant preserves each comparison outcome for every active lane entering the helper.',
             'Bounds checks remain in the shared opcode template, now using a uniform compile-time register limit. Original generic fallback and outer structural validation are byte-identical. No dependence on merely presumed external admission is introduced.',
             'This removes varying metadata from the instruction helper ABI; whether it materially improves execution is unproven. At most the existing opcode set per distinct register count is instantiated, with no per-instruction template explosion.',
             'Compile or execute only if Z fails its decisive original frame and its terminal evidence has been preserved. If Z succeeds, discard this unused contingent source.'
         ]}
publish(W / 'ispc-derivation.json', (json.dumps(proof, indent=2) + '\n').encode())
print(json.dumps({'source_stages': 3, 'instructions': instruction_count, 'layers': layer_count,
                  'compiled': False, 'executed_workload': False}))
