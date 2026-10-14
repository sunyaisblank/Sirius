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

    def origin(value):
        seen = set()
        while words(value)[0] == "OpCopyObject":
            assert value not in seen
            seen.add(value)
            value = words(value)[2]
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

    same = lambda value: lambda other: origin(value) == origin(other)
    number = lambda n: lambda value: constant(value, n)
    mask = lambda value: binary(value, "OpBitwiseAnd", same(hi), number(0xfff00000), True)
    range_condition = one(key for key in local if binary(key, "OpIEqual", mask, number(0x41900000), True))
    fraction = lambda value: binary(value, "OpBitwiseAnd", same(hi), number(0xfffff), True)
    hidden = lambda value: binary(value, "OpBitwiseOr", fraction, number(0x100000), True)
    upper = lambda value: binary(value, "OpShiftLeftLogical", hidden, number(6))
    lower = lambda value: binary(value, "OpShiftRightLogical", same(lo), number(26))
    guess = one(key for key in local if binary(key, "OpBitwiseOr", upper, lower, True))
    product = one(key for key in local if words(key)[0] == "OpFunctionCall"
                  and names.get(words(key)[2]) == "PB32WideProduct"
                  and [origin(x) for x in words(key)[3:]] == [origin(guess), d])
    plo = one(key for key in local if words(key)[:3] == ["OpCompositeExtract", u32, product] and words(key)[3:] == ["0"])
    phi = one(key for key in local if words(key)[:3] == ["OpCompositeExtract", u32, product] and words(key)[3:] == ["1"])
    nlo = one(key for key in local if binary(key, "OpShiftLeftLogical", same(m), number(26)))
    nhi = one(key for key in local if binary(key, "OpShiftRightLogical", same(m), number(6)))
    below = lambda value: isop(value, "OpULessThan", phi, nhi)
    equal = lambda value: binary(value, "OpIEqual", same(phi), same(nhi), True)
    low_ok = lambda value: isop(value, "OpULessThanEqual", plo, nlo)
    equal_low_ok = lambda value: binary(value, "OpLogicalAnd", equal, low_ok, True)
    not_above = one(key for key in local if binary(key, "OpLogicalOr", below, equal_low_ok, True))
    rlo = one(key for key in local if isop(key, "OpISub", nlo, plo))
    high_difference = lambda value: isop(value, "OpISub", nhi, phi)
    borrow_test = one(key for key in local if isop(key, "OpULessThan", nlo, plo))
    borrow = lambda value: (words(origin(value))[0] == "OpSelect" and
                            origin(words(origin(value))[2]) == origin(borrow_test) and
                            constant(words(origin(value))[3], 1) and constant(words(origin(value))[4], 0))
    rhi = one(key for key in local if binary(key, "OpISub", high_difference, borrow))
    high_zero = lambda value: binary(value, "OpIEqual", same(rhi), number(0), True)
    low_small = lambda value: isop(value, "OpULessThan", rlo, d)
    verified = one(key for key in local if binary(key, "OpLogicalAnd", high_zero, low_small, True))
    nonzero = lambda value: binary(value, "OpINotEqual", same(rlo), number(0), True)
    sticky = lambda value: (words(origin(value))[0] == "OpSelect" and
                            nonzero(words(origin(value))[2]) and
                            constant(words(origin(value))[3], 1) and constant(words(origin(value))[4], 0))
    jam = one(key for key in local if binary(key, "OpBitwiseOr", same(guess), sticky, True))
    rounded = one(key for key in local if words(key)[0] == "OpFunctionCall"
                  and names.get(words(key)[2]) == "PB32Round" and origin(words(key)[-1]) == origin(jam))

    # Exact nested single-entry true arms. No guess/product/remainder can bypass
    # its preceding acceptance predicate, and success returns this round call.
    blocks, location, current = {}, {}, None
    for line in body.splitlines():
        line = line.strip()
        match = re.match(r"(%\S+) = OpLabel$", line)
        if match:
            current = match[1]
            blocks[current] = []
        if current is not None:
            blocks[current].append(line)
            match = re.match(r"(%\S+) =", line)
            if match:
                location[match[1]] = current
    predecessors = {label: set() for label in blocks}
    for label, lines in blocks.items():
        for line in lines:
            tokens = line.split()
            targets = tokens[1:] if tokens[:1] == ["OpBranch"] else tokens[2:] if tokens[:1] == ["OpBranchConditional"] else []
            for target in targets:
                predecessors[target].add(label)
    def true_arm(condition, expected):
        label = location[condition]
        branches = [line.split() for line in blocks[label] if line.startswith("OpBranchConditional ")]
        branch = one(branches)
        assert branch[1] == condition and branch[2] == expected and predecessors[expected] == {label}
    true_arm(range_condition, location[product])
    assert location[guess] == location[product]
    true_arm(not_above, location[rlo])
    true_arm(verified, location[rounded])
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
            "guess": guess, "wide_product_call": product, "exact_remainder": [rlo, rhi],
            "acceptance_predicates": [range_condition, not_above, verified],
            "integer_range_before_guess_extraction": True,
            "exact_product_remainder_inequalities": True, "exact_remainder_sticky": True,
            "accepted_path_returns_same_round_sign_exponent": True,
            "original_complete_fallback_suffix_normalized_equal": True,
            "unsigned64_or_float_to_uint_required": False,
            "arithmetic_equivalence_executed": False, "cost_or_hit_rate_established": False}
