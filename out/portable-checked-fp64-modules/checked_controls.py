def portable_normal_sum_controls(text):
    """Keep original binary32 checks; isolate the three untrusted seed bodies."""
    if sorted(re.findall(r'OpCapability (\S+)',text)) != ['Float64','Shader']:
        raise ValueError('checked estimates require precisely Shader/Float64')
    if set(re.findall(r'OpTypeInt (\d+) [01]',text)) != {'32'}:
        raise ValueError('checked estimates must keep integer arithmetic32-only')
    wide, = re.findall(r'(%\S+) = OpTypeFloat 64',text)
    filtered=text
    bodies={}
    for name in ['PB32Divide','PB32Sqrt','PB32SeedUnsigned']:
        found=re.findall(r'%'+name+r'(?:_\d+)? = OpFunction .*?OpFunctionEnd',text,re.S)
        if len(found)!=1:raise ValueError('seed body not uniquely retained: '+name)
        bodies[name]=found[0];filtered=filtered.replace(found[0],'')
    divide,sqrt,seed=[bodies[n] for n in ['PB32Divide','PB32Sqrt','PB32SeedUnsigned']]
    if re.findall(r'= (OpFMul|OpFDiv) '+re.escape(wide)+r' ',divide)!=['OpFMul','OpFDiv']:
        raise ValueError('division estimate changed its two native operations')
    if re.findall(r'= (OpFMul|OpFDiv) '+re.escape(wide)+r' ',sqrt)!=['OpFMul'] or len(re.findall(r'= OpExtInst '+re.escape(wide)+r' %\S+ Sqrt ',sqrt))!=1:
        raise ValueError('square-root estimate changed its two native operations')
    if divide.count('OpLoopMerge')!=1 or sqrt.count('OpLoopMerge')!=1:
        raise ValueError('unchanged integer fallbacks must remain present')
    if seed.count('OpConvertFToU')!=1 or 'OpPhi '+wide not in seed or 'OpBitcast %v2uint' not in seed or 'OpCompositeExtract %uint' not in seed:
        raise ValueError('seed conversion lost its guarded raw-bit sanitization')
    if re.findall(r'= OpULessThan %bool %\S+ (%\S+)',seed)!=['%uint_1102053376'] or not re.search(r'%uint_1102053376 = OpConstant %uint 1102053376',text):
        raise ValueError('seed conversion lost positive finite <2^28 bound')
    filtered=re.sub(r'^.*OpCapability Float64\n','',filtered,flags=re.M)
    filtered=re.sub(r'^.*= OpTypeFloat 64\n','',filtered,flags=re.M)
    # This call checks every original native32 operation, dependency chain,
    # NoContraction, conversion exclusion and float mode without modification.
    stock_portable_normal_sum_controls(filtered)
    entry=re.search(r'OpEntryPoint GLCompute (%\S+)',text)[1]
    text=text.replace('OpCapability Shader','OpCapability Shader\n               OpCapability RoundingModeRTE',1)
    local=re.search(r'^.*OpExecutionMode '+re.escape(entry)+r' LocalSize.*$',text,re.M)[0]
    return text.replace(local,local+'\n               OpExecutionMode '+entry+' RoundingModeRTE 32',1)
