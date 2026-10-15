"""Whole actual-module source controls; only the future owned reader runs these."""
import re


def run(stock,candidate,function_body,inspect):
    body=function_body(candidate,'PB32Divide')
    def line(old,new):
        pattern=r'(?m)^(\s*)'+re.escape(old)+r'$'
        assert len(re.findall(pattern,body))==1,old
        changed=re.sub(pattern,lambda m:m[1]+new,body)
        assert changed!=body and candidate.count(body)==1
        return candidate.replace(body,changed,1)
    def pair(first,second):
        text=line(*first);old,new=second
        pattern=r'(?m)^(\s*)'+re.escape(old)+r'$'
        assert len(re.findall(pattern,text))==1,old
        return re.sub(pattern,lambda m:m[1]+new,text)
    negatives=[
      ('coordinated_fast_and_fallback_wrong_exponent','typed_CFG_semantic',pair(
        ('%3019 = OpFunctionCall %uint %PB32Round %sign_2 %2943 %3018','%3019 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %3018'),
        ('%3047 = OpFunctionCall %uint %PB32Round %sign_2 %2943 %3046','%3047 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %3046'))),
      ('selector_phi_arms_swapped','typed_CFG_semantic',line('%2948 = OpPhi %uint %uint_1291845632 %2946 %uint_1283457024 %2947','%2948 = OpPhi %uint %uint_1283457024 %2946 %uint_1291845632 %2947')),
      ('selector_other_dominating_predicate','typed_CFG_semantic',line('OpBranchConditional %2939 %2946 %2947','OpBranchConditional %2892 %2946 %2947')),
      ('proposal_operands_swapped','typed_CFG_semantic',line('%proposal = OpFDiv %float %2951 %2952','%proposal = OpFDiv %float %2952 %2951')),
      ('numerator_uses_denominator_bits','typed_CFG_semantic',line('%2951 = OpBitcast %float %numeratorBits','%2951 = OpBitcast %float %denominatorBits')),
      ('range_mask_wrong','typed_CFG_semantic',line('%2954 = OpBitwiseAnd %uint %2953 %uint_4286578688','%2954 = OpBitwiseAnd %uint %2953 %uint_8388607')),
      ('range_binade_wrong','typed_CFG_semantic',line('%2955 = OpIEqual %bool %2954 %uint_1283457024','%2955 = OpIEqual %bool %2954 %uint_1291845632')),
      ('h_shift_wrong','typed_CFG_semantic',line('%h = OpShiftLeftLogical %uint %2959 %int_3','%h = OpShiftLeftLogical %uint %2959 %int_1')),
      ('guess_margin_wrong','typed_CFG_semantic',line('%guess = OpISub %uint %h %uint_16','%guess = OpISub %uint %h %uint_8')),
      ('wide_product_wrong_divisor','typed_CFG_semantic',line('%2960 = OpUMulExtended %PB32Wide %guess %2938','%2960 = OpUMulExtended %PB32Wide %guess %remainder_1')),
      ('wide_high_uses_low','typed_CFG_semantic',line('%2962 = OpCompositeExtract %uint %2960 1','%2962 = OpCompositeExtract %uint %2960 0')),
      ('borrow_inclusive','typed_CFG_semantic',line('%2964 = OpULessThan %bool %numeratorLo %2961','%2964 = OpULessThanEqual %bool %numeratorLo %2961')),
      ('borrow_phi_arms_swapped','typed_CFG_semantic',line('%2968 = OpPhi %uint %uint_1 %2966 %uint_0 %2967','%2968 = OpPhi %uint %uint_0 %2966 %uint_1 %2967')),
      ('borrow_added','typed_CFG_semantic',line('%residualHi = OpISub %uint %2963 %2968','%residualHi = OpIAdd %uint %2963 %2968')),
      ('high_zero_inverted','typed_CFG_semantic',line('%2969 = OpIEqual %bool %residualHi %uint_0','%2969 = OpINotEqual %bool %residualHi %uint_0')),
      ('residual_bound_wrong','typed_CFG_semantic',line('%2973 = OpShiftLeftLogical %uint %2938 %int_5','%2973 = OpShiftLeftLogical %uint %2938 %int_3')),
      ('residual_bound_inclusive','typed_CFG_semantic',line('%2974 = OpULessThan %bool %residualLo %2973','%2974 = OpULessThanEqual %bool %residualLo %2973')),
      ('residual_guard_false_becomes_true','typed_CFG_semantic',line('%2975 = OpPhi %bool %2974 %2971 %false %2972','%2975 = OpPhi %bool %2974 %2971 %true %2972')),
      ('restoring_branch_bypasses_guard','typed_CFG_semantic',line('OpBranchConditional %2975 %2977 %2976','OpBranchConditional %2955 %2977 %2976')),
      ('restoring_final_quotient_subtracts','typed_CFG_semantic',line('%exactQuotient = OpIAdd %uint %guess %digit_3','%exactQuotient = OpISub %uint %guess %digit_3')),
      ('restoring_quotient_uses_partial_digit','typed_CFG_semantic',line('%exactQuotient = OpIAdd %uint %guess %digit_3','%exactQuotient = OpIAdd %uint %guess %digit_2')),
      ('sticky_uses_pre_correction_residual','typed_CFG_semantic',line('%3013 = OpINotEqual %bool %exactRemainder_3 %uint_0','%3013 = OpINotEqual %bool %residualLo %uint_0')),
      ('sticky_uses_borrow_predicate','typed_CFG_semantic',line('OpBranchConditional %3013 %3015 %3016','OpBranchConditional %2964 %3015 %3016')),
      ('sticky_phi_arms_swapped','typed_CFG_semantic',line('%3017 = OpPhi %uint %uint_1 %3015 %uint_0 %3016','%3017 = OpPhi %uint %uint_0 %3015 %uint_1 %3016')),
      ('original_special_guard_operand_changed','typed_CFG_semantic',line('%2892 = OpFunctionCall %bool %PB32Nan %a_32','%2892 = OpFunctionCall %bool %PB32Nan %b_17')),
      ('original_fallback_loop_bound_changed','typed_CFG_semantic',line('%3032 = OpULessThan %bool %3026 %uint_4','%3032 = OpULessThan %bool %3026 %uint_3')),
      ('original_fallback_sticky_lost','typed_CFG_semantic',line('%quotient = OpBitwiseOr %uint %3021 %uint_1','%quotient = OpBitwiseOr %uint %3021 %uint_0')),
      ('original_fallback_divisor_changed','typed_CFG_semantic',line('%digit_4 = OpUDiv %uint %numerator_0 %2938','%digit_4 = OpUDiv %uint %numerator_0 %remainder_1')),
      ('original_fallback_round_exponent_changed','typed_CFG_semantic',line('%3047 = OpFunctionCall %uint %PB32Round %sign_2 %2943 %3046','%3047 = OpFunctionCall %uint %PB32Round %sign_2 %exponent_7 %3046')),
      ('phi_predecessor_mismatch','parser_closed_CFG',line('%2968 = OpPhi %uint %uint_1 %2966 %uint_0 %2967','%2968 = OpPhi %uint %uint_1 %2966 %uint_0 %2947')),
      ('switch_side_entrance','parser_closed_CFG',line('OpSwitch %int_0 %3031','OpSwitch %int_0 %3031 1 %2977')),
      ('extra_typed_call_in_selector_arm','typed_CFG_semantic',line('%2946 = OpLabel','%2946 = OpLabel\n%graph_extra = OpFunctionCall %bool %PB32Zero %a_32')),
      ('extra_typed_arithmetic_in_restoring_arm','typed_CFG_semantic',line('%threshold4 = OpShiftLeftLogical %uint %2938 %int_4','%threshold4 = OpShiftLeftLogical %uint %2938 %int_4\n%graph_extra = OpIAdd %uint %residualLo %uint_1')),
      ('extra_corrective_UDiv_in_restoring_arm','typed_CFG_semantic',line('%threshold4 = OpShiftLeftLogical %uint %2938 %int_4','%threshold4 = OpShiftLeftLogical %uint %2938 %int_4\n%graph_extra = OpUDiv %uint %residualLo %2938')),
    ]
    # Each row names actual SSA from a distinct restoring position. All mutated
    # operands remain uint/bool and dominate their uses; no tool validation claimed.
    rows=[
      (4,'%residualLo','%threshold4','%2978','%2982','%2980','%2981','%2983','%exactRemainder','%2984','%digit','%uint_0'),
      (3,'%exactRemainder','%threshold3','%2985','%2989','%2987','%2988','%2990','%exactRemainder_0','%2991','%digit_0','%digit'),
      (2,'%exactRemainder_0','%threshold2','%2992','%2996','%2994','%2995','%2997','%exactRemainder_1','%2998','%digit_1','%digit_0'),
      (1,'%exactRemainder_1','%threshold1','%2999','%3003','%3001','%3002','%3004','%exactRemainder_2','%3005','%digit_2','%digit_1'),
      (0,'%exactRemainder_2','%2938','%3006','%3010','%3008','%3009','%3011','%exactRemainder_3','%3012','%digit_3','%digit_2'),
    ]
    for s,prior,threshold,test,mask,yes,no,selected,remainder,bit,digit,old_digit in rows:
        prefix='restoring_step'+str(s)+'_'
        negatives.extend([
          (prefix+'inclusive_direction_wrong','typed_CFG_semantic',line(test+' = OpUGreaterThanEqual %bool '+prior+' '+threshold,test+' = OpULessThan %bool '+prior+' '+threshold)),
          (prefix+'mask_arms_swapped','typed_CFG_semantic',line(mask+' = OpPhi %uint %uint_4294967295 '+yes+' %uint_0 '+no,mask+' = OpPhi %uint %uint_0 '+yes+' %uint_4294967295 '+no)),
          (prefix+'subtract_adds','typed_CFG_semantic',line(remainder+' = OpISub %uint '+prior+' '+selected,remainder+' = OpIAdd %uint '+prior+' '+selected)),
          (prefix+'wrong_selected_mask','typed_CFG_semantic',line(selected+' = OpBitwiseAnd %uint '+threshold+' '+mask,selected+' = OpBitwiseAnd %uint '+threshold+' %2968')),
          (prefix+'digit_bit_wrong','typed_CFG_semantic',line(bit+' = OpBitwiseAnd %uint %uint_'+str(1<<s)+' '+mask,bit+' = OpBitwiseAnd %uint %uint_'+str(1<<(s-1) if s else 1+1)+' '+mask)),
        ])
        if s:
            negatives.append((prefix+'threshold_wrong_shift','typed_CFG_semantic',line(threshold+' = OpShiftLeftLogical %uint %2938 %int_'+str(s),threshold+' = OpShiftLeftLogical %uint %2938 %int_'+str(s-1))))
        else:
            negatives.append((prefix+'threshold_wrong_divisor','typed_CFG_semantic',line(test+' = OpUGreaterThanEqual %bool '+prior+' %2938',test+' = OpUGreaterThanEqual %bool '+prior+' %remainder_1')))
        if s<4:
            negatives.append((prefix+'uses_stale_residual','typed_CFG_semantic',line(test+' = OpUGreaterThanEqual %bool '+prior+' '+threshold,test+' = OpUGreaterThanEqual %bool %residualLo '+threshold)))
            negatives.append((prefix+'drops_previous_digit','typed_CFG_semantic',line(digit+' = OpBitwiseOr %uint '+old_digit+' '+bit,digit+' = OpBitwiseOr %uint %uint_0 '+bit)))
    annotation='               OpDecorate %proposal NoContraction'
    assert candidate.count(annotation)==1
    negatives.extend([
      ('unsigned_NoWrap_on_modular_subtraction','annotation_contract',candidate.replace(annotation,'               OpDecorate %residualLo NoUnsignedWrap\n'+annotation,1)),
      ('signed_NoWrap_on_modular_low_subtraction','annotation_contract',candidate.replace(annotation,'               OpDecorate %residualLo NoSignedWrap\n'+annotation,1)),
      ('proposal_NoContraction_missing','annotation_contract',candidate.replace(annotation+'\n','',1)),
    ])
    assert all(text!=candidate for _,_,text in negatives)
    positives=[('exact_actual_module',candidate)]
    ids=re.findall(r'^\s*(%\S+) =',body,re.M);aliases={value:'%graph_renamed_'+str(i) for i,value in enumerate(ids)}
    positives.append(('all_Divide_local_IDs_renamed_in_full_module',re.sub(r'%[^\s]+',lambda m:aliases.get(m[0],m[0]),candidate)))
    positives.append(('equivalent_leading_whitespace_and_blank_lines','\n\n'+'\n'.join('  '+ln.strip() if ln.strip() else '' for ln in candidate.splitlines())+'\n\n'))
    accepted=[]
    for label,text in positives:
        proof=inspect(stock,text,function_body)
        assert proof['original_complete_prefix_normalized_equal'] and proof['original_complete_fallback_suffix_normalized_equal']
        assert proof['corrective_UDiv_absent'] and len(proof['restoring_steps'])==5
        accepted.append(label)
    rejected=[]
    for label,scope,text in negatives:
        try:inspect(stock,text,function_body)
        except (AssertionError,KeyError,ValueError,IndexError):rejected.append({'name':label,'scope':scope})
        else:raise AssertionError('required negative accepted: '+label)
    return {'schema':'sirius-fp32-estimate-restoring-fullmodule-graph-controls-v1','passed':True,
            'positive_count':len(accepted),'negative_count':len(rejected),'accepted':accepted,'rejected':rejected,
            'mutants_are_complete_actual_module':True,'all_five_restoring_positions_mutated':True,
            'typed_CFG_semantic_mutations':sum(scope=='typed_CFG_semantic' for _,scope,_ in negatives),
            'mutated_modules_assembled_or_validated':False,'scope':'Source graph-reader controls only; operand-type/CFG mutations are not independently SPIR-V-validated. Parser-closed-CFG controls are explicitly separate. No execution, numerical or cost admission.'}
