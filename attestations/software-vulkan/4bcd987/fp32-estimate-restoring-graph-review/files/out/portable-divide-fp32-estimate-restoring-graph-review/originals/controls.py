"""Whole actual-module source controls; run only by the future owned reader."""
import re


def run(stock,candidate,function_body,inspect):
    body=function_body(candidate,'PB32Divide')
    def line(old,new):
        # Full instruction matching avoids the historical uint_21/uint_214... prefix collision.
        pattern=r'(?m)^(\s*)'+re.escape(old)+r'$'
        assert len(re.findall(pattern,body))==1,old
        changed=re.sub(pattern,lambda m:m[1]+new,body)
        assert changed!=body and candidate.count(body)==1
        return candidate.replace(body,changed,1)
    coordinated=line('%2972 = OpFunctionCall %uint %PB32Round %sign_2 %2930 %2971','%2972 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %2971')
    old_tail='%3000 = OpFunctionCall %uint %PB32Round %sign_2 %2930 %2999'
    pattern=r'(?m)^(\s*)'+re.escape(old_tail)+r'$'
    assert len(re.findall(pattern,coordinated))==1
    coordinated=re.sub(pattern,lambda m:m[1]+'%3000 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %2999',coordinated)
    negatives=[
      ('coordinated_fast_and_fallback_wrong_exponent','typed_CFG_semantic',coordinated),
      ('selector_phi_arms_swapped','typed_CFG_semantic',line('%2935 = OpPhi %uint %uint_1291845632 %2933 %uint_1283457024 %2934','%2935 = OpPhi %uint %uint_1283457024 %2933 %uint_1291845632 %2934')),
      ('selector_uses_other_dominating_bool','typed_CFG_semantic',line('OpBranchConditional %2926 %2933 %2934','OpBranchConditional %2879 %2933 %2934')),
      ('proposal_operands_swapped','typed_CFG_semantic',line('%proposal = OpFDiv %float %2938 %2939','%proposal = OpFDiv %float %2939 %2938')),
      ('numerator_bitcast_uses_denominator','typed_CFG_semantic',line('%2938 = OpBitcast %float %numeratorBits','%2938 = OpBitcast %float %denominatorBits')),
      ('range_mask_wrong','typed_CFG_semantic',line('%2941 = OpBitwiseAnd %uint %2940 %uint_4286578688','%2941 = OpBitwiseAnd %uint %2940 %uint_8388607')),
      ('range_binade_wrong','typed_CFG_semantic',line('%2942 = OpIEqual %bool %2941 %uint_1283457024','%2942 = OpIEqual %bool %2941 %uint_1291845632')),
      ('h_shift_wrong','typed_CFG_semantic',line('%h = OpShiftLeftLogical %uint %2946 %int_3','%h = OpShiftLeftLogical %uint %2946 %int_1')),
      ('guess_margin_wrong','typed_CFG_semantic',line('%guess = OpISub %uint %h %uint_16','%guess = OpISub %uint %h %uint_8')),
      ('wide_product_wrong_divisor','typed_CFG_semantic',line('%2947 = OpUMulExtended %PB32Wide %guess %2925','%2947 = OpUMulExtended %PB32Wide %guess %remainder_1')),
      ('wide_product_high_uses_low','typed_CFG_semantic',line('%2949 = OpCompositeExtract %uint %2947 1','%2949 = OpCompositeExtract %uint %2947 0')),
      ('borrow_inclusive_comparison','typed_CFG_semantic',line('%2951 = OpULessThan %bool %numeratorLo %2948','%2951 = OpULessThanEqual %bool %numeratorLo %2948')),
      ('borrow_phi_arms_swapped','typed_CFG_semantic',line('%2955 = OpPhi %uint %uint_1 %2953 %uint_0 %2954','%2955 = OpPhi %uint %uint_0 %2953 %uint_1 %2954')),
      ('borrow_added_instead_of_subtracted','typed_CFG_semantic',line('%residualHi = OpISub %uint %2950 %2955','%residualHi = OpIAdd %uint %2950 %2955')),
      ('high_zero_predicate_inverted','typed_CFG_semantic',line('%2956 = OpIEqual %bool %residualHi %uint_0','%2956 = OpINotEqual %bool %residualHi %uint_0')),
      ('residual_bound_wrong_shift','typed_CFG_semantic',line('%2960 = OpShiftLeftLogical %uint %2925 %int_5','%2960 = OpShiftLeftLogical %uint %2925 %int_3')),
      ('residual_bound_inclusive','typed_CFG_semantic',line('%2961 = OpULessThan %bool %residualLo %2960','%2961 = OpULessThanEqual %bool %residualLo %2960')),
      ('residual_guard_phi_false_becomes_true','typed_CFG_semantic',line('%2962 = OpPhi %bool %2961 %2958 %false %2959','%2962 = OpPhi %bool %2961 %2958 %true %2959')),
      ('correction_branch_bypasses_exact_guard','typed_CFG_semantic',line('OpBranchConditional %2962 %2964 %2963','OpBranchConditional %2942 %2964 %2963')),
      ('correction_wrong_numerator','typed_CFG_semantic',line('%digit = OpUDiv %uint %residualLo %2925','%digit = OpUDiv %uint %numeratorLo %2925')),
      ('correction_wrong_divisor','typed_CFG_semantic',line('%digit = OpUDiv %uint %residualLo %2925','%digit = OpUDiv %uint %residualLo %remainder_1')),
      ('correction_low_product_wrong_digit','typed_CFG_semantic',line('%2965 = OpIMul %uint %digit %2925','%2965 = OpIMul %uint %guess %2925')),
      ('exact_quotient_subtracts_digit','typed_CFG_semantic',line('%exactQuotient = OpIAdd %uint %guess %digit','%exactQuotient = OpISub %uint %guess %digit')),
      ('sticky_uses_borrow_predicate','typed_CFG_semantic',line('OpBranchConditional %2966 %2968 %2969','OpBranchConditional %2951 %2968 %2969')),
      ('sticky_phi_arms_swapped','typed_CFG_semantic',line('%2970 = OpPhi %uint %uint_1 %2968 %uint_0 %2969','%2970 = OpPhi %uint %uint_0 %2968 %uint_1 %2969')),
      ('round_uses_pre_normalized_exponent','typed_CFG_semantic',line('%2972 = OpFunctionCall %uint %PB32Round %sign_2 %2930 %2971','%2972 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %2971')),
      ('original_special_guard_operand_changed','typed_CFG_semantic',line('%2879 = OpFunctionCall %bool %PB32Nan %a_32','%2879 = OpFunctionCall %bool %PB32Nan %b_17')),
      ('original_fallback_loop_bound_changed','typed_CFG_semantic',line('%2985 = OpULessThan %bool %2979 %uint_4','%2985 = OpULessThan %bool %2979 %uint_3')),
      ('original_fallback_sticky_lost','typed_CFG_semantic',line('%quotient = OpBitwiseOr %uint %2974 %uint_1','%quotient = OpBitwiseOr %uint %2974 %uint_0')),
      ('original_fallback_divisor_changed','typed_CFG_semantic',line('%digit_0 = OpUDiv %uint %numerator_0 %2925','%digit_0 = OpUDiv %uint %numerator_0 %remainder_1')),
      ('original_fallback_round_exponent_changed','typed_CFG_semantic',line('%3000 = OpFunctionCall %uint %PB32Round %sign_2 %2930 %2999','%3000 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %2999')),
      ('phi_predecessor_mismatch','parser_closed_CFG',line('%2955 = OpPhi %uint %uint_1 %2953 %uint_0 %2954','%2955 = OpPhi %uint %uint_1 %2953 %uint_0 %2934')),
      ('switch_side_entrance','parser_closed_CFG',line('OpSwitch %int_0 %2984','OpSwitch %int_0 %2984 1 %2964')),
      ('extra_typed_call_in_selector_arm','typed_CFG_semantic',line('%2933 = OpLabel','%2933 = OpLabel\n%graph_extra = OpFunctionCall %bool %PB32Zero %a_32')),
      ('extra_typed_arithmetic_in_fast_block','typed_CFG_semantic',line('%h = OpShiftLeftLogical %uint %2946 %int_3','%h = OpShiftLeftLogical %uint %2946 %int_3\n%graph_extra = OpIAdd %uint %h %uint_1')),
      ('unsigned_NoWrap_on_modular_subtraction','annotation_contract',candidate.replace('               OpDecorate %proposal NoContraction','               OpDecorate %residualLo NoUnsignedWrap\n               OpDecorate %proposal NoContraction',1)),
      ('signed_NoWrap_on_modular_low_subtraction','annotation_contract',candidate.replace('               OpDecorate %proposal NoContraction','               OpDecorate %residualLo NoSignedWrap\n               OpDecorate %proposal NoContraction',1)),
      ('proposal_NoContraction_missing','annotation_contract',candidate.replace('               OpDecorate %proposal NoContraction\n','',1)),
    ]
    assert all(text!=candidate for _,_,text in negatives)
    positives=[('exact_actual_module',candidate)]
    ids=re.findall(r'^\s*(%\S+) =',body,re.M)
    aliases={value:'%graph_renamed_'+str(i) for i,value in enumerate(ids)}
    renamed=re.sub(r'%[^\s]+',lambda m:aliases.get(m[0],m[0]),candidate)
    positives.append(('all_Divide_local_IDs_renamed_in_full_module',renamed))
    formatted='\n\n'+'\n'.join('  '+ln.strip() if ln.strip() else '' for ln in candidate.splitlines())+'\n\n'
    positives.append(('equivalent_leading_whitespace_and_blank_lines',formatted))
    accepted=[]
    for label,text in positives:
        proof=inspect(stock,text,function_body)
        assert proof['original_complete_prefix_normalized_equal'] and proof['original_complete_fallback_suffix_normalized_equal']
        accepted.append(label)
    rejected=[]
    for label,scope,text in negatives:
        try:inspect(stock,text,function_body)
        except (AssertionError,KeyError,ValueError,IndexError):rejected.append({'name':label,'scope':scope})
        else:raise AssertionError('required negative accepted: '+label)
    return {'schema':'sirius-fp32-estimate-fullmodule-graph-controls-v1','passed':True,
            'positive_count':len(accepted),'negative_count':len(rejected),'accepted':accepted,'rejected':rejected,
            'mutants_are_complete_actual_module':True,'typed_CFG_semantic_mutations':sum(scope=='typed_CFG_semantic' for _,scope,_ in negatives),
            'mutated_modules_assembled_or_validated':False,'scope':'Source graph-reader rejection controls; typed/CFG mutations preserve declared operand types and edges unless explicitly parser_closed_CFG. No execution or numerical admission.'}
