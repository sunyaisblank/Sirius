"""Accepted product reader with bounded O1 camera extension; no build/runtime entry."""
import re

def function_body(text, name):
    match = re.search(r'OpName (%\S+) "' + re.escape(name) + '"', text)
    assert match, name
    found = re.search(r"^\s*" + re.escape(match[1]) + r" = OpFunction\b.*?OpFunctionEnd",
                      text, re.M | re.S)
    assert found, name
    return found.group(0) + "\n"


def normalized_function(text, name):
    body = function_body(text, name)
    definitions = dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))
    names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', text))
    local = {identifier: "local" + str(i) for i, identifier in enumerate(re.findall(r"^\s*(%\S+) =", body, re.M))}
    def value(identifier):
        if identifier in local:
            return local[identifier]
        expression = definitions[identifier]
        if expression.startswith("OpFunction "):
            return ["function", names[identifier]]
        assert expression.startswith(("OpType", "OpConstant", "OpSpecConstant")), expression
        return [value(word) if word.startswith("%") else word for word in expression.split()]
    return [[value(word) if word.startswith("%") else word for word in line.split()]
            for line in body.splitlines() if line.strip()]


def function_slot_choice(body, definitions, value):
    """Recognize only a private uint load after one exact two-store diamond."""
    load = definitions[value].split()
    assert len(load) == 3 and load[0] == "OpLoad"
    kind, slot = load[1:]
    assert definitions[kind] == "OpTypeInt 32 0"
    variable = definitions[slot].split()
    assert len(variable) == 3 and variable[0] == "OpVariable" and variable[2] == "Function"
    assert definitions[variable[1]].split() == ["OpTypePointer", "Function", kind]
    blocks, rows, current = {}, [], None
    for line in body.splitlines():
        words = line.split()
        if not words:
            continue
        result, op, args = ((words[0], words[2], words[3:]) if len(words) > 1 and words[1] == "="
                            else (None, words[0], words[1:]))
        row = {"result": result, "op": op, "args": args}
        rows.append(row)
        if op == "OpLabel":
            current = result
            assert current not in blocks
            blocks[current] = []
        if current is not None and op != "OpFunctionEnd":
            blocks[current].append(row)
    # The slot has no aliases, address arithmetic, pointer copies, calls or stores as data.
    for row in rows:
        if slot in row["args"]:
            assert ((row["op"] == "OpStore" and len(row["args"]) == 2 and row["args"][0] == slot and row["args"][1] != slot) or
                    (row["op"] == "OpLoad" and row["result"] is not None and row["args"] == [kind, slot]))
    predecessors = {name: set() for name in blocks}
    for name, block in blocks.items():
        tail = block[-1]
        if tail["op"] == "OpBranch":
            assert len(tail["args"]) == 1
            targets = tail["args"]
        elif tail["op"] == "OpBranchConditional":
            assert len(tail["args"]) == 3
            targets = tail["args"][1:]
        else:
            assert tail["op"] in ("OpReturnValue", "OpReturn", "OpUnreachable")
            targets = []
        for target in targets:
            predecessors[target].add(name)
    positions = [(name, index) for name, block in blocks.items()
                 for index, row in enumerate(block) if row["result"] == value]
    assert len(positions) == 1
    merge, index = positions[0]
    # Each arm ends immediately after its store. Nothing in the merge may rewrite the slot before this load.
    assert all(not (row["op"] == "OpStore" and row["args"][0] == slot)
               for row in blocks[merge][:index])
    candidates = []
    for header, block in blocks.items():
        if len(block) < 3 or block[-2]["op"] != "OpSelectionMerge" or block[-2]["args"] != [merge, "None"]:
            continue
        branch = block[-1]
        assert branch["op"] == "OpBranchConditional" and len(branch["args"]) == 3
        condition, true_block, false_block = branch["args"]
        assert true_block != false_block
        assert predecessors[true_block] == predecessors[false_block] == {header}
        assert predecessors[merge] == {true_block, false_block}
        values = []
        for arm in (true_block, false_block):
            exact = blocks[arm]
            assert len(exact) == 3 and exact[0]["op"] == "OpLabel"
            assert exact[1]["op"] == "OpStore" and len(exact[1]["args"]) == 2 and exact[1]["args"][0] == slot
            assert exact[2]["op"] == "OpBranch" and exact[2]["args"] == [merge]
            values.append(exact[1]["args"][1])
        candidates.append({"slot": slot, "header": header, "merge": merge, "condition": condition,
                           "true_block": true_block, "false_block": false_block,
                           "yes": values[0], "no": values[1], "load": value,
                           "direct_private_slot_only": True, "no_side_entrances_or_intervening_store": True})
    assert len(candidates) == 1
    return candidates[0]


def multiply_shape(text, camera=False, allow_o0=False):
    """Require exact full product/shift/jam graph; Camera O1 has8inlined generic sites."""
    body = (re.search(r"^\s*" + re.escape(re.search(r"OpEntryPoint GLCompute (%\S+)", text)[1]) + r" = OpFunction\b.*?OpFunctionEnd", text, re.M | re.S)[0] if camera else function_body(text, "PB32Multiply"))
    definitions = dict(re.findall(r"^\s*(%\S+) = (.+)$", text, re.M))
    names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', text))
    integers = dict(re.findall(r"(%\S+) = OpTypeInt (\d+) [01]", text))
    wide = [identifier for identifier, width in integers.items() if width == "64"]
    assert len(wide) == 1 and definitions[wide[0]] == "OpTypeInt 64 0"
    wide = wide[0]
    words = lambda identifier: definitions[identifier].split()
    def origin(identifier):
        seen = set()
        while words(identifier)[0] == "OpCopyObject":
            assert identifier not in seen
            seen.add(identifier); identifier = words(identifier)[2]
        return identifier
    def number(identifier):
        identifier = origin(identifier); expression = words(identifier)
        if expression[0] == "OpConstant":
            return int(expression[2], 0)
        if expression[0] in ("OpShiftLeftLogical", "OpISub"):
            a, b = number(expression[2]), number(expression[3])
            return a << b if expression[0] == "OpShiftLeftLogical" else a - b
        raise ValueError("not an accepted integer constant expression: " + identifier)
    def constant(identifier, expected):
        try: return number(identifier) == expected
        except ValueError: return False
    loads = {}
    def load_choice(identifier):
        identifier = origin(identifier)
        assert allow_o0 and not camera
        if identifier not in loads:
            loads[identifier] = function_slot_choice(body, definitions, identifier)
        return loads[identifier]
    def same_value(left, right):
        left, right = origin(left), origin(right)
        if left == right:
            return True
        if not allow_o0 or words(left)[0] != "OpLoad" or words(right)[0] != "OpLoad":
            return False
        epoch = lambda item: {key: value for key, value in item.items() if key != "load"}
        return epoch(load_choice(left)) == epoch(load_choice(right))
    def choice(identifier, condition, yes, no):
        expression = words(origin(identifier))
        if expression[0] == "OpSelect":
            return origin(expression[2]) == condition and constant(expression[3], yes) and constant(expression[4], no)
        if expression[0] == "OpLoad" and allow_o0:
            proof = load_choice(identifier)
            return (origin(proof["condition"]) == condition and
                    constant(proof["yes"], yes) and constant(proof["no"], no))
        if expression[0] != "OpPhi" or len(expression) != 6:
            return False
        incoming = {expression[i + 1]: expression[i] for i in (2, 4)}
        branches = re.findall(r"OpBranchConditional " + re.escape(condition) + r" (%\S+) (%\S+)", body)
        return any(set(incoming) == {a, b} and constant(incoming[a], yes) and constant(incoming[b], no) for a, b in branches)
    products = re.findall(r"(%\S+) = OpIMul " + re.escape(wide) + r" (%\S+) (%\S+)", body)
    expected = 8 if camera else 1
    assert len(products) == expected
    assert body.count("OpUMulExtended") == (1 if camera else 0)
    assert len(re.findall(r"= OpIMul " + re.escape(wide) + r" ", text)) == expected
    results = []
    for product, left, right in products:
        for operand in (left, right):
            expression = words(origin(operand)); assert expression[:2] == ["OpUConvert", wide]
            extraction = words(origin(expression[2]))
            if camera:
                assert integers.get(extraction[1]) == "32"
            else:
                assert extraction[0] == "OpCompositeExtract" and extraction[-1] == "0"
                unpack = words(origin(extraction[2])); assert unpack[0] == "OpFunctionCall" and names[unpack[2]] == "PB32Unpack"
        top = [identifier for identifier, expression in definitions.items()
               if expression.startswith("OpBitwiseAnd " + wide + " ") and
               ((origin(words(identifier)[2]) == product and constant(words(identifier)[3], 1 << 47)) or
                (origin(words(identifier)[3]) == product and constant(words(identifier)[2], 1 << 47)))]
        assert len(top) == 1
        def nonzero_of(value):
            results = []
            for identifier, expression in definitions.items():
                parts = expression.split()
                if parts[0] == "OpINotEqual" and len(parts) == 4 and (
                        (origin(parts[2]) == value and constant(parts[3], 0)) or
                        (origin(parts[3]) == value and constant(parts[2], 0))): results.append(identifier)
            assert len(results) == 1
            return results[0]
        upper = nonzero_of(top[0])
        shifts = re.findall(r"(%\S+) = OpShiftRightLogical " + re.escape(wide) + " " + re.escape(product) + r" (%\S+)", body)
        assert len(shifts) == 1
        shifted, distance = shifts[0]
        assert choice(distance, upper, 21, 20)
        masks = []
        for identifier, expression in definitions.items():
            parts = expression.split()
            if parts[:2] != ["OpISub", wide] or not constant(parts[3], 1): continue
            shift = words(origin(parts[2]))
            if shift[:2] == ["OpShiftLeftLogical", wide] and constant(shift[2], 1) and same_value(shift[3], distance): masks.append(identifier)
        assert len(masks) == 1
        discarded = [identifier for identifier, expression in definitions.items() if
                     expression.startswith("OpBitwiseAnd " + wide + " ") and
                     {origin(word) for word in words(identifier)[2:]} == {product, masks[0]}]
        assert len(discarded) == 1
        sticky_test = nonzero_of(discarded[0])
        narrow = [identifier for identifier, expression in definitions.items() if
                  expression.startswith("OpUConvert ") and origin(words(identifier)[2]) == shifted and integers[words(identifier)[1]] == "32"]
        assert len(narrow) == 1
        extensions = []
        for identifier in re.findall(r"^\s*(%\S+) = OpBitwiseOr .*", body, re.M):
            operands = [origin(word) for word in words(identifier)[2:]]
            if narrow[0] in operands and choice(operands[1] if operands[0] == narrow[0] else operands[0], sticky_test, 1, 0): extensions.append(identifier)
        assert len(extensions) == 1
        if not camera:
            rounded = [identifier for identifier in re.findall(r"^\s*(%\S+) = OpFunctionCall .*", body, re.M)
                       if names.get(words(identifier)[2]) == "PB32Round" and origin(words(identifier)[-1]) == extensions[0]]
            assert len(rounded) == 1 and re.search(r"OpReturnValue " + re.escape(rounded[0]) + r"\b", body)
        else:
            # Original Round source is frozen unchanged; O1 inlines it. Do not invent an out-of-line call witness.
            assert len(re.findall(re.escape(extensions[0]) + r"(?![A-Za-z0-9_])", body)) > 1
        results.append({"wide_product": product, "upper_bit": 47, "shift_values": [20, 21],
                        "complete_discarded_bit_jam": True, "narrow32_sticky_or": extensions[0]})
    return {"unsigned64_product_sites": expected, "inlined_camera": camera,
            "normalized_unsigned_operands": not camera, "uint32_operands_widened": True,
            "complete_product_shift_jam_graphs": results,
            "function_slot_load_choices": list(loads.values()),
            "common_round_scope": "Original source unchanged; out-of-line call joined except O1 camera where Round remains inlined."}
