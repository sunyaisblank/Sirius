"""Inspect the actual named production functions and op5-specialized CFG."""
from pathlib import Path
import hashlib
import json
import re

work = Path(__file__).resolve().parent
build = Path.cwd() / 'bin/linux-gcc/portable-denominator-preparation'


def parse(text):
    constants = {}
    functions = {}
    current = None
    for line in text.splitlines():
        match = re.match(r'\s*(?:(%\S+)\s*=\s*)?(Op\w+)\s*(.*)', line)
        if not match:
            continue
        result, opcode, rest = match.groups()
        args = rest.split()
        if opcode == 'OpConstant' and len(args) == 2:
            try:
                constants[result] = int(args[1])
            except ValueError:
                pass
        elif opcode in ('OpConstantTrue', 'OpConstantFalse'):
            constants[result] = opcode == 'OpConstantTrue'
        if opcode == 'OpFunction':
            current = []
            functions[result] = current
        if current is not None:
            current.append((result, opcode, args))
        if opcode == 'OpFunctionEnd':
            current = None
    return constants, functions


def blocks(instructions):
    result = {}
    current = None
    for item in instructions:
        if item[1] == 'OpLabel':
            current = item[0]
            result[current] = []
        elif current is not None:
            result[current].append(item)
    return result


def successors(items, known):
    for _, opcode, args in reversed(items):
        if opcode == 'OpBranch':
            return [args[0]]
        if opcode == 'OpBranchConditional':
            value = known.get(args[0])
            return args[1:3] if value is None else [args[1] if value else args[2]]
        if opcode == 'OpSwitch':
            value = known.get(args[0])
            if value is None:
                return [args[1], *args[3::2]]
            choices = {int(args[i]): args[i + 1] for i in range(2, len(args), 2)}
            return [choices.get(value, args[1])]
        if opcode in ('OpReturn', 'OpReturnValue', 'OpUnreachable'):
            return []
    raise AssertionError('Missing terminator')


def reachable(entry, graph):
    seen = set()
    pending = [entry]
    while pending:
        node = pending.pop()
        if node not in seen:
            seen.add(node)
            pending.extend(graph[node])
    return seen


def specialize_op5(instructions, constants):
    cfg = blocks(instructions)
    entry = next(iter(cfg))
    parameters = [r for r, op, _ in instructions if op == 'OpFunctionParameter']
    assert len(parameters) == 3
    known = dict(constants)
    known[parameters[2]] = 5
    seen = set(cfg)
    for _ in range(len(cfg) + 1):
        previous = (set(seen), dict(known))
        for block, items in cfg.items():
            if block not in seen:
                continue
            for result, opcode, args in items:
                value = None
                if opcode in ('OpIEqual', 'OpINotEqual'):
                    a, b = (known.get(x) for x in args[1:3])
                    if a is not None and b is not None:
                        value = (a == b) if opcode == 'OpIEqual' else (a != b)
                elif opcode == 'OpLogicalNot':
                    a = known.get(args[1])
                    if a is not None:
                        value = not a
                elif opcode == 'OpPhi':
                    values = [known.get(args[i]) for i in range(1, len(args), 2)
                              if args[i + 1] in seen]
                    if values and all(v is not None and v == values[0] for v in values):
                        value = values[0]
                if value is not None:
                    assert result not in known or known[result] == value
                    known[result] = value
        graph = {node: successors(items, known) for node, items in cfg.items()}
        seen = reachable(entry, graph)
        if previous == (seen, known):
            break
    else:
        raise AssertionError('CFG specialization did not converge')
    graph = {n: [s for s in graph[n] if s in seen] for n in seen}
    dominators = {n: ({entry} if n == entry else set(seen)) for n in seen}
    for _ in range(len(seen) + 1):
        old = {n: set(v) for n, v in dominators.items()}
        for node in seen - {entry}:
            parents = [p for p, edges in graph.items() if node in edges]
            assert parents
            dominators[node] = {node} | set.intersection(*(dominators[p] for p in parents))
        if old == dominators:
            break
    else:
        raise AssertionError('Dominators did not converge')
    return cfg, graph, seen, dominators


def calls(items, function):
    return [item for item in items if item[1] == 'OpFunctionCall' and item[2][1] == function]


def call_closure(start, functions):
    seen = set()
    pending = [start]
    while pending:
        f = pending.pop()
        if f not in seen:
            seen.add(f)
            pending.extend(args[1] for _, op, args in functions[f] if op == 'OpFunctionCall')
    return seen


records = []
for stem in ('retained_transport_portable', 'retained_transport_portable_normal_sum',
             'retained_endpoint_portable_normal_sum', 'retained_dopri_phase_portable_normal_sum'):
    path = build / (stem + '.spvasm')
    text = path.read_text()
    constants, functions = parse(text)
    body = functions['%RTEvaluateArithmetic']
    cfg, graph, seen, dom = specialize_op5(body, constants)
    prep = calls(body, '%RPPrepareDivideSeedDenominator')
    seeds = calls(body, '%RPDivideSeedPrepared')
    assert len(prep) == 1 and len(seeds) == 2
    locate = lambda item: next(node for node, items in cfg.items() if item in items)
    prep_block = locate(prep[0])
    seed_blocks = [locate(seed) for seed in seeds]
    assert prep_block in seen and all(prep_block in dom[node] for node in seed_blocks)
    assert all(prep_block not in reachable(s, graph) for s in graph[prep_block])
    loops = [node for node, items in cfg.items() if node in seen
             and any(op == 'OpLoopMerge' for _, op, _ in items)]
    assert len(loops) == 1 and prep_block in dom[loops[0]]
    assert all(not calls(cfg[node], '%RPDivideSeed') for node in seen)
    assert any(calls(items, '%RPDivideSeed') for node, items in cfg.items() if node not in seen)
    stores = [args for _, op, args in body if op == 'OpStore' and args[0] == '%prepared']
    assert stores == [['%prepared', prep[0][0]]]
    initial_pointer = seeds[0][2][3]
    assert any(op == 'OpStore' and args == [initial_pointer, prep[0][0]] for _, op, args in body)
    correction_pointer = seeds[1][2][3]
    loaded = [r for r, op, args in body if op == 'OpLoad' and args[1] == '%prepared']
    assert len(loaded) == 1
    assert any(op == 'OpStore' and args == [correction_pointer, loaded[0]] for _, op, args in body)
    prepare = functions['%RPPrepareDivideSeedDenominator']
    divisions = [args for _, op, args in prepare if op == 'OpUDiv']
    assert len(divisions) == 1 and constants[divisions[0][1]] == 0xffffffff
    closure = call_closure('%RPDivideSeedPrepared', functions)
    assert '%PB32DivideNormalized' in closure and '%PB32Divide' not in closure
    assert all(op != 'OpUDiv' for f in closure for _, op, _ in functions[f])
    seed_body = functions['%RPDivideSeedPrepared']
    scaling = calls(seed_body, '%OrderedMultiply32')
    powers = calls(seed_body, '%RPPowerOfTwo')
    assert len(scaling) == 2 and len(powers) == 1
    definitions = {r: (op, args) for r, op, args in seed_body if r}
    min_normal_op, min_normal_args = definitions[scaling[0][2][3]]
    assert min_normal_op == 'OpFunctionCall' and min_normal_args[1] == '%RPFromBits'
    assert constants[min_normal_args[2]] == 0x00800000
    scaled_op, scaled_args = definitions[scaling[1][2][2]]
    assert scaled_op == 'OpPhi'
    assert set(scaled_args[1::2]) == {scaling[0][0], scaling[0][2][2]}
    assert scaling[1][2][3] == powers[0][0]
    exponent_op, exponent_args = definitions[powers[0][2][2]]
    assert exponent_op == 'OpPhi'
    additions = [(r, args) for r, op, args in seed_body if op == 'OpIAdd'
                 and constants.get(args[2]) == 126]
    assert len(additions) == 1
    assert set(exponent_args[1::2]) == {additions[0][0], additions[0][1][1]}
    normalized = functions['%PB32DivideNormalized']
    assert sum(op == 'OpUMulExtended' for _, op, _ in normalized) == 1
    assert len(calls(normalized, '%PB32Round')) == 1
    assert sum(op == 'OpLoopMerge' for _, op, _ in normalized) == 1
    assert any(op == 'OpUGreaterThanEqual' for _, op, _ in normalized)
    module = build / 'products' / (stem + '.spv')
    records.append({'module': stem, 'module_sha256': hashlib.sha256(module.read_bytes()).hexdigest(),
                    'assembly_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                    'op5_reachable_blocks': len(seen), 'phase_loop': loops[0],
                    'single_preparation_block': prep_block, 'prepared_seed_blocks': seed_blocks,
                    'preparation_dominates_op5_phase_loop_and_both_seed_sites': True,
                    'preparation_outside_cycles': True, 'single_owned_value_reused': True,
                    'no_generic_seed_on_specialized_op5_route': True,
                    'generic_seed_retained_on_sqrt_route': True,
                    'prepared_seed_call_closure_has_no_UDiv': True,
                    'both_original_scaling_operations_retained': True})
report = {'passed': True, 'records': records,
          'scope': 'Exact generated production function CFG with op=5 branches specialized. '
                   'Source argument establishes phases/correction count and numerical equivalence; '
                   'this finite compiler screen is not execution, timing or machine-checked universal proof.'}
(work / 'compiled-route-check.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'passed': True, 'actual_production_modules': len(records)}))
