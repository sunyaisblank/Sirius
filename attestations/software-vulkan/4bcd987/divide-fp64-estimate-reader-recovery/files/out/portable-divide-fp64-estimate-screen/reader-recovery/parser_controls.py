"""Bounded full-module text negatives; never assemble, compile or execute SPIR-V."""


def controls(stock, candidate, graph, function_body):
    def replace_line(text, before, after):
        lines = text.splitlines(keepends=True)
        selected = [i for i, line in enumerate(lines) if line.strip() == before]
        assert len(selected) == 1, before
        index = selected[0]
        indent = lines[index][:len(lines[index]) - len(lines[index].lstrip())]
        lines[index] = "".join(indent + row + "\n" for row in after.splitlines())
        return "".join(lines)

    def change(before, after):
        return lambda text: replace_line(text, before, after)

    cases = [
        ("range_mask", change("%2939 = OpBitwiseAnd %uint %2938 %uint_4293918720", "%2939 = OpBitwiseAnd %uint %2938 %uint_2147483647")),
        ("range_exponent", change("%2940 = OpIEqual %bool %2939 %uint_1099956224", "%2940 = OpIEqual %bool %2939 %uint_0")),
        ("swapped_member_store", change("OpStore %2934 %2937", "OpStore %2934 %2938")),
        ("second_member_store", change("OpStore %2934 %2937", "OpStore %2934 %2937\nOpStore %2934 %2937")),
        ("pointer_alias", change("%2935 = OpAccessChain %_ptr_Function_uint %estimateBits %int_1", "%2935 = OpAccessChain %_ptr_Function_uint %estimateBits %int_1\n%reader_bad_alias = OpCopyObject %_ptr_Function_uint %2934")),
        ("load_result_type", change("%2943 = OpLoad %uint %2935", "%2943 = OpLoad %int %2935")),
        ("switch_side_entrance", change("OpSwitch %int_0 %2998", "OpSwitch %int_0 %2998 0 %2942")),
        ("range_true_side_entrance", change("%2931 = OpLabel", "%2931 = OpLabel\nOpBranch %2942")),
        ("reversed_product_bound_predicate", change("OpBranchConditional %2951 %2953 %2954", "OpBranchConditional %2951 %2954 %2953")),
        ("wrong_product_low_bound", change("%2960 = OpULessThanEqual %bool %2959 %numeratorLo", "%2960 = OpULessThanEqual %bool %2959 %numeratorHi")),
        ("wrong_or_phi_arm", change("%2962 = OpPhi %bool %true %2953 %2961 %2956", "%2962 = OpPhi %bool %false %2953 %2961 %2956")),
        ("wrong_product_high_limb", change("%2950 = OpCompositeExtract %uint %2949 1", "%2950 = OpCompositeExtract %uint %2949 0")),
        ("wrong_product_denominator", change("%2949 = OpUMulExtended %PB32Wide %guess %2927", "%2949 = OpUMulExtended %PB32Wide %guess %2933")),
        ("wrong_numerator_high_shift", change("%numeratorHi = OpShiftRightLogical %uint %2933 %int_6", "%numeratorHi = OpShiftRightLogical %uint %2933 %int_26")),
        ("wrong_borrow_phi", change("%2971 = OpPhi %uint %uint_1 %2969 %uint_0 %2970", "%2971 = OpPhi %uint %uint_0 %2969 %uint_1 %2970")),
        ("omitted_borrow", change("%residualHi = OpISub %uint %2966 %2971", "%residualHi = OpIAdd %uint %2966 %2971")),
        ("wrong_remainder_bound", change("%2976 = OpULessThan %bool %residualLo %2927", "%2976 = OpULessThan %bool %residualLo %numeratorLo")),
        ("wrong_and_phi_arm", change("%2977 = OpPhi %bool %2976 %2974 %false %2975", "%2977 = OpPhi %bool %2976 %2974 %true %2975")),
        ("wrong_sticky_phi", change("%2984 = OpPhi %uint %uint_1 %2982 %uint_0 %2983", "%2984 = OpPhi %uint %uint_0 %2982 %uint_1 %2983")),
        ("wrong_round_exponent", change("%2986 = OpFunctionCall %uint %PB32Round %sign_2 %2932 %2985", "%2986 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %2985")),
    ]

    def load_before_store(text):
        text = replace_line(text, "%2943 = OpLoad %uint %2935", "")
        return replace_line(text, "OpStore %2935 %2938", "%2943 = OpLoad %uint %2935\nOpStore %2935 %2938")

    def nested_side_entrance(text):
        return replace_line(text, "%2958 = OpLabel", "%2958 = OpLabel\nOpBranch %2952")

    def fallback_digit(text):
        rows = [line.strip() for line in function_body(text, "PB32Divide").splitlines() if " = OpUDiv " in line]
        assert len(rows) == 1
        tokens = rows[0].split()
        return replace_line(text, rows[0], " ".join(tokens[:-1] + ["%uint_1"]))

    cases.extend([("load_before_unique_store", load_before_store),
                  ("nested_phi_side_entrance", nested_side_entrance),
                  ("changed_original_fallback", fallback_digit)])
    positive = graph.inspect(stock, candidate, function_body)
    rows = [{"name": "exact_complete_original_candidate", "expected_acceptance": True,
             "actual_acceptance": True}]
    for name, mutation in cases:
        modified = mutation(candidate)
        assert modified != candidate
        try:
            graph.inspect(stock, modified, function_body)
        except (AssertionError, KeyError, ValueError, IndexError):
            accepted = False
        else:
            accepted = True
        assert not accepted, name
        rows.append({"name": name, "expected_acceptance": False, "actual_acceptance": accepted})
    return {"passed": True, "controls": rows, "negative_count": len(cases),
            "positive_graph": positive,
            "scope": "Complete immutable candidate disassembly text only. Some pointer/CFG negatives are deliberately malformed parser inputs; no shader validity, execution or numerical claim. Original module tools/validation are reused unchanged."}
