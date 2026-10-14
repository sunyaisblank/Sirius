"""Strict one-function SSA screen; no arithmetic or provider execution."""
import re


def inspect(stock, candidate, function_body):
    def module(text):
        return dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))

    definitions = module(candidate)
    names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', candidate))
    body = function_body(candidate, "PB32Divide")
    local = module(body)
    assert len(local) < 400

    def words(value):
        return definitions[value].split()

    # Function-local CFG only, bounded by this complete Divide body.
    blocks, location, positions, current = {}, {}, {}, None
    for line in body.splitlines():
        line = line.strip()
        match = re.match(r"(%\S+) = OpLabel$", line)
        if match:
            current = match[1]; blocks[current] = []
        if current is not None:
            positions[(current, line)] = len(blocks[current]); blocks[current].append(line)
            match = re.match(r"(%\S+) =", line)
            if match: location[match[1]] = current
    predecessors = {label: set() for label in blocks}
    for label, lines in blocks.items():
        for line in lines:
            tokens = line.split()
            targets = (tokens[1:] if tokens[:1] == ["OpBranch"] else
                       tokens[2:4] if tokens[:1] == ["OpBranchConditional"] else
                       [tokens[2], *tokens[4::2]] if tokens[:1] == ["OpSwitch"] else [])
            for target in targets: predecessors[target].add(label)
    aliases, choices, slot_proof = {}, {}, {}

    def origin(value):
        seen = set()
        while value in aliases or words(value)[0] == "OpCopyObject":
            assert value not in seen
            seen.add(value)
            value = aliases[value] if value in aliases else words(value)[2]
        return value

    def isop(value, opcode, *operands):
        row = words(origin(value))
        return row[0] == opcode and [origin(x) for x in row[2:]] == [origin(x) for x in operands]

    def constant(value, expected):
        row = words(origin(value))
        return row[0] == "OpConstant" and float(row[2]) == expected

    def one(values):
        values = list(values)
        assert len(values) == 1, values
        return values[0]

    wide = one(key for key, value in definitions.items() if value == "OpTypeFloat 64")
    u32 = one(key for key, value in definitions.items() if value == "OpTypeInt 32 0")
    vector = one(key for key, value in definitions.items() if value == "OpTypeVector " + u32 + " 2")
    divides = [key for key in local if words(key)[:2] == ["OpFDiv", wide]]
    divide = one(divides)
    scaled, denominator = words(divide)[2:]
    assert words(origin(scaled))[:2] == ["OpFMul", wide]
    conversion, scale = words(origin(scaled))[2:]
    assert constant(scale, 67108864)
    assert words(origin(conversion))[:2] == ["OpConvertUToF", wide]
    assert words(origin(denominator))[:2] == ["OpConvertUToF", wide]
    m, d = origin(words(origin(conversion))[2]), origin(words(origin(denominator))[2])
    bits = one(key for key in local if isop(key, "OpBitcast", divide) and words(key)[1] == vector)
    lo = one(key for key in local if words(key)[:3] == ["OpCompositeExtract", u32, bits] and words(key)[3:] == ["0"])
    hi = one(key for key in local if words(key)[:3] == ["OpCompositeExtract", u32, bits] and words(key)[3:] == ["1"])

    def binary(value, opcode, predicate_a, predicate_b, commute=False):
        row = words(origin(value))
        if row[0] != opcode or len(row) != 4:
            return False
        a, b = row[2:]
        return (predicate_a(a) and predicate_b(b)) or (commute and predicate_a(b) and predicate_b(a))

    def choice(value, condition, yes, no):
        value = origin(value); row = words(value)
        if row[0] == "OpSelect":
            if condition(row[2]) and yes(row[3]) and no(row[4]):
                choices[value] = {"entry": location[value], "merge": location[value], "kind": "Select"}
                return True
            return False
        if row[0] != "OpPhi" or len(row) != 6:
            return False
        merge = location[value]
        headers = [label for label, lines in blocks.items() if "OpSelectionMerge " + merge + " None" in lines]
        if len(headers) != 1: return False
        header = headers[0]
        branches = [line.split() for line in blocks[header] if line.startswith("OpBranchConditional ")]
        if len(branches) != 1 or len(branches[0]) != 4: return False
        branch = branches[0]
        if not condition(branch[1]): return False
        incoming = {row[i + 1]: row[i] for i in (2, 4)}
        if len(incoming) != 2 or predecessors[merge] != set(incoming): return False
        def arm(destination, predicate):
            matches = [label for label, incoming_value in incoming.items() if predicate(incoming_value)]
            if len(matches) != 1: return False
            label = matches[0]; incoming_value = origin(incoming[label])
            if label == destination:
                return predecessors[destination] == {header} and blocks[destination][-1] == "OpBranch " + merge
            nested = choices.get(incoming_value)
            return (nested is not None and nested["entry"] == destination and nested["merge"] == label
                    and predecessors[destination] == {header} and blocks[label][-1] == "OpBranch " + merge)
        if not arm(branch[2], yes) or not arm(branch[3], no): return False
        choices[value] = {"entry": header, "merge": merge, "kind": "closed two-arm Phi", "predicate": branch[1]}
        return True

    def boolean_constant(value, expected):
        return words(origin(value))[0] == ("OpConstantTrue" if expected else "OpConstantFalse")

    def logical(value, opcode, left, right):
        if binary(value, opcode, left, right, True):
            choices[origin(value)] = {"entry": location[origin(value)], "merge": location[origin(value)], "kind": opcode}
            return True
        if opcode == "OpLogicalAnd":
            return choice(value, left, right, lambda x: boolean_constant(x, False))
        assert opcode == "OpLogicalOr"
        return choice(value, left, lambda x: boolean_constant(x, True), right)

    same = lambda value: lambda other: origin(value) == origin(other)
    number = lambda n: lambda value: constant(value, n)
    mask = lambda value: binary(value, "OpBitwiseAnd", same(hi), number(0xfff00000), True)
    range_condition = one(key for key in local if binary(key, "OpIEqual", mask, number(0x41900000), True))
    # Only the two uniquely written private uint members may resolve a load.
    # The stores precede the range branch; every load is in its single-entry
    # true block. All Function-pointer uses are direct Store/Load; no escape,
    # alias, additional store, side entry or global substitution is admitted.
    range_header = location[range_condition]
    range_branch = [line.split() for line in blocks[range_header] if line.startswith("OpBranchConditional ")]
    assert len(range_branch) == 1 and range_branch[0][1] == range_condition
    range_true = range_branch[0][2]
    loads = [key for key in local if words(key)[0] == "OpLoad"]
    if loads:
        pointer_rows = [words(words(key)[2]) for key in loads]
        assert all(row[0] == "OpAccessChain" and len(row) == 4 for row in pointer_rows)
        roots = {row[2] for row in pointer_rows}; assert len(roots) == 1
        variable = roots.pop(); variable_row = words(variable)
        assert variable_row[0] == "OpVariable" and variable_row[2:] == ["Function"]
        pointer_type = words(variable_row[1]); assert pointer_type[:2] == ["OpTypePointer", "Function"]
        assert words(pointer_type[2]) == ["OpTypeStruct", u32, u32]
        members = [key for key in local if words(key)[0] == "OpAccessChain" and words(key)[2] == variable]
        assert len(members) == 2
        bindings = {}
        all_stores = [line for line in body.splitlines() if line.strip().startswith("OpStore ")]
        assert len(all_stores) == 2
        for member in members:
            row = words(member); assert words(row[1]) == ["OpTypePointer", "Function", u32]
            index = 0 if constant(row[3], 0) else 1 if constant(row[3], 1) else None
            assert index is not None and index not in bindings
            stores = [line.strip() for line in all_stores if line.split()[1] == member]
            assert len(stores) == 1 and origin(stores[0].split()[2]) == origin(lo if index == 0 else hi)
            assert stores[0] in blocks[range_header] and location[member] == range_header
            condition_line = next(line for line in blocks[range_header] if line.startswith(range_condition + " ="))
            assert positions[(range_header, stores[0])] < positions[(range_header, condition_line)]
            for line in body.splitlines():
                tokens = line.split()
                if member not in tokens: continue
                opcode = tokens[2] if len(tokens) > 2 and tokens[1] == "=" else tokens[0]
                assert tokens[0] == member or opcode in ("OpStore", "OpLoad"), line
            bindings[index] = member
        for line in body.splitlines():
            tokens = line.split()
            if variable not in tokens: continue
            assert tokens[0] == variable or (len(tokens) > 2 and tokens[2] == "OpAccessChain"), line
        for value in loads:
            assert words(value)[:2] == ["OpLoad", u32]
            member = words(value)[2]
            assert member in bindings.values() and location[value] == range_true and predecessors[range_true] == {range_header}
            index = next(i for i, pointer in bindings.items() if pointer == member)
            aliases[value] = lo if index == 0 else hi
        slot_proof = {"variable": variable, "member_pointers": bindings,
                      "load_ids": loads, "single_stores_dominate_single_entry_guarded_loads": True,
                      "direct_only_no_alias_escape_or_extra_store": True}
    else:
        assert not re.search(r"\bOpStore\b", body)
    fraction = lambda value: binary(value, "OpBitwiseAnd", same(hi), number(0xfffff), True)
    hidden = lambda value: binary(value, "OpBitwiseOr", fraction, number(0x100000), True)
    upper = lambda value: binary(value, "OpShiftLeftLogical", hidden, number(6))
    lower = lambda value: binary(value, "OpShiftRightLogical", same(lo), number(26))
    guess = one(key for key in local if binary(key, "OpBitwiseOr", upper, lower, True))
    def exact_product(key):
        row = words(key)
        if row[0] == "OpUMulExtended": return [origin(x) for x in row[2:]] == [origin(guess), d]
        return (row[0] == "OpFunctionCall" and names.get(row[2]) == "PB32WideProduct"
                and [origin(x) for x in row[3:]] == [origin(guess), d])
    product = one(key for key in local if exact_product(key))
    assert words(words(product)[1]) == ["OpTypeStruct", u32, u32]
    def member(index):
        values = [key for key in local if words(key)[:3] == ["OpCompositeExtract", u32, product]
                  and words(key)[3:] == [str(index)]]
        assert 1 <= len(values) <= 2
        for value in values[1:]: aliases[value] = values[0]
        return values[0]
    plo, phi = member(0), member(1)
    nlo = one(key for key in local if binary(key, "OpShiftLeftLogical", same(m), number(26)))
    nhi = one(key for key in local if binary(key, "OpShiftRightLogical", same(m), number(6)))
    below = lambda value: isop(value, "OpULessThan", phi, nhi)
    equal = lambda value: binary(value, "OpIEqual", same(phi), same(nhi), True)
    low_ok = lambda value: isop(value, "OpULessThanEqual", plo, nlo)
    equal_low_ok = lambda value: logical(value, "OpLogicalAnd", equal, low_ok)
    not_above = one(key for key in local if logical(key, "OpLogicalOr", below, equal_low_ok))
    rlo = one(key for key in local if isop(key, "OpISub", nlo, plo))
    high_difference = lambda value: isop(value, "OpISub", nhi, phi)
    borrow_test = one(key for key in local if isop(key, "OpULessThan", nlo, plo))
    borrow = lambda value: choice(value, same(borrow_test), number(1), number(0))
    rhi = one(key for key in local if binary(key, "OpISub", high_difference, borrow))
    high_zero = lambda value: binary(value, "OpIEqual", same(rhi), number(0), True)
    low_small = lambda value: isop(value, "OpULessThan", rlo, d)
    verified = one(key for key in local if logical(key, "OpLogicalAnd", high_zero, low_small))
    nonzero = lambda value: binary(value, "OpINotEqual", same(rlo), number(0), True)
    sticky = lambda value: choice(value, nonzero, number(1), number(0))
    jam = one(key for key in local if binary(key, "OpBitwiseOr", same(guess), sticky, True))
    rounded = one(key for key in local if words(key)[0] == "OpFunctionCall"
                  and names.get(words(key)[2]) == "PB32Round" and origin(words(key)[-1]) == origin(jam))

    # Exact closed control paths from range gate through product, borrow,
    # remainder gate and sticky. No side entrance can bypass an acceptance gate.
    def true_arm(condition, expected):
        label = location[condition]
        branches = [line.split() for line in blocks[label] if line.startswith("OpBranchConditional ")]
        branch = one(branches)
        assert branch[1] == condition and branch[2] == expected and predecessors[expected] == {label}
    true_arm(range_condition, location[product])
    assert location[guess] == location[product]
    assert choices[origin(not_above)]["entry"] == location[product]
    true_arm(not_above, location[rlo])
    borrow_value = origin(words(rhi)[3])
    assert choices[borrow_value]["entry"] == location[rlo]
    assert choices[borrow_value]["merge"] == location[rhi]
    assert choices[origin(verified)]["entry"] == location[rhi]
    jam_words = words(jam)
    sticky_value = origin(jam_words[3] if origin(jam_words[2]) == origin(guess) else jam_words[2])
    true_arm(verified, choices[sticky_value]["entry"])
    assert choices[sticky_value]["merge"] == location[rounded]
    assert blocks[location[rounded]][-1] == "OpReturnValue " + rounded

    # Compare the complete original loop/sticky/round/return suffix after ID
    # normalization. Only its single preheader label and unchanged m/d/sign/e
    # inputs are abstracted; constants, widths, edges and every instruction stay.
    def fallback(text):
        all_defs = module(text)
        local_names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', text))
        fragment = function_body(text, "PB32Divide")
        defs = module(fragment)
        ids = lambda key: all_defs[key].split()
        udiv = one(key for key in defs if ids(key)[0] == "OpUDiv")
        numerator, divisor = ids(udiv)[2:]
        assert ids(numerator)[0] == "OpShiftLeftLogical"
        remainder = ids(numerator)[2]
        assert ids(remainder)[0] == "OpPhi"
        loop = one(re.finditer(r"OpLoopMerge (%\S+) (%\S+) (\S+)", fragment))
        labels = list(re.finditer(r"^\s*(%\S+) = OpLabel$", fragment[:loop.start()], re.M))
        start = labels[-1].start()
        suffix = fragment[start:]
        local_ids = re.findall(r"^\s*(%\S+) =", suffix, re.M)
        mapped = {value: "local" + str(i) for i, value in enumerate(local_ids)}
        incoming = ids(remainder)[2:]
        initial = one(incoming[i] for i in range(0, len(incoming), 2) if incoming[i + 1] not in mapped)
        assert ids(initial)[:1] == ["OpISub"] and ids(initial)[3] == divisor
        normalized_m = ids(initial)[2]
        calls = [ids(key) for key in defs if ids(key)[0] == "OpFunctionCall"
                 and local_names.get(ids(key)[2]) == "PB32Round"]
        tail_call = one(row for row in calls if row[-1] in mapped)
        roles = {normalized_m: "normalized_m", divisor: "divisor", tail_call[3]: "sign", tail_call[4]: "exponent"}
        def value(identifier):
            if identifier in mapped:
                return mapped[identifier]
            if identifier in roles:
                return roles[identifier]
            row = ids(identifier)
            if row[0] == "OpLabel":
                return "single_preheader"
            if row[0] == "OpFunction":
                return ["function", dict(re.findall(r'OpName (%\S+) "([^"]*)"', text))[identifier]]
            assert row[0].startswith(("OpType", "OpConstant")) or row[0] == "OpISub", row
            return [value(word) if word.startswith("%") else word for word in row]
        normalized = [[value(word) if word.startswith("%") else word for word in line.split()]
                      for line in suffix.splitlines() if line.strip()]
        return normalized, tail_call, normalized_m, divisor
    before, _, _, _ = fallback(stock)
    after, tail_call, fallback_m, fallback_d = fallback(candidate)
    assert before == after and origin(fallback_m) == m and origin(fallback_d) == d
    assert [origin(x) for x in words(rounded)[3:5]] == [origin(x) for x in tail_call[3:5]]
    return {"native_estimate_accuracy_assumed": False, "wide_division": divide,
            "guess": guess, "exact_product": {"id": product, "opcode": words(product)[0]}, "exact_remainder": [rlo, rhi],
            "acceptance_predicates": [range_condition, not_above, verified],
            "integer_range_before_guess_extraction": True,
            "exact_product_remainder_inequalities": True, "exact_remainder_sticky": True,
            "accepted_path_returns_same_round_sign_exponent": True,
            "original_complete_fallback_suffix_normalized_equal": True,
            "unsigned64_or_float_to_uint_required": False,
            "Function_uint2_member_load_proof": slot_proof, "closed_choice_proofs": choices,
            "named_PB32WideProduct_required": False,
            "arithmetic_equivalence_executed": False, "cost_or_hit_rate_established": False}
