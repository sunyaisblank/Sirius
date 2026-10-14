"""Bounded inlined-camera estimate footprint; frozen-source Round is explicit."""
import re


def inspect(stock, candidate):
    definitions = dict(re.findall(r'^\s*(%\S+) = (.+)$', candidate, re.M))
    word = lambda value: definitions[value].split()
    def one(values):
        values = list(values)
        assert len(values) == 1, values
        return values[0]
    u32 = one(k for k,v in definitions.items() if v == 'OpTypeInt 32 0')
    wide = one(k for k,v in definitions.items() if v == 'OpTypeFloat 64')
    vec = one(k for k,v in definitions.items() if v == 'OpTypeVector '+u32+' 2')
    entry = re.search(r'OpEntryPoint GLCompute (%\S+)', candidate)[1]
    body = re.search(r'^\s*'+re.escape(entry)+r' = OpFunction\b.*?OpFunctionEnd', candidate, re.M|re.S)[0]
    local = dict(re.findall(r'^\s*(%\S+) = (.+)$',body,re.M))
    assert not re.search(r'OpName %\S+ "PB32Divide"',candidate)
    assert not re.search(r'OpName %\S+ "PB32Round"',candidate)
    blocks,where,current = {},{},None
    for line in body.splitlines():
        line=line.strip(); match=re.fullmatch(r'(%\S+) = OpLabel',line)
        if match: current=match[1];blocks[current]=[]
        if current is not None:
            blocks[current].append(line); match=re.match(r'(%\S+) =',line)
            if match:where[match[1]]=current
    pred={k:set() for k in blocks}
    for label,lines in blocks.items():
        for line in lines:
            t=line.split(); targets=(t[1:] if t[:1]==['OpBranch'] else t[2:4] if t[:1]==['OpBranchConditional'] else [t[2],*t[4::2]] if t[:1]==['OpSwitch'] else [])
            for target in targets:pred[target].add(label)
    aliases={}
    def origin(value):
        seen=set()
        while value in aliases or word(value)[0]=='OpCopyObject':
            assert value not in seen;seen.add(value);value=aliases[value] if value in aliases else word(value)[2]
        return value
    def number(value,n):
        t=word(origin(value));return t[0]=='OpConstant' and float(t[2])==n
    def binary(value,op,a,b,commute=False):
        t=word(origin(value))
        return t[0]==op and len(t)==4 and ((a(t[2]) and b(t[3])) or (commute and a(t[3]) and b(t[2])))
    same=lambda x:lambda y:origin(x)==origin(y)
    num=lambda n:lambda x:number(x,n)
    choices={}
    def choice(value,condition,yes,no):
        value=origin(value);t=word(value)
        if t[0]=='OpSelect':
            if condition(t[2]) and yes(t[3]) and no(t[4]):
                choices[value]={'entry':where[value],'merge':where[value]};return True
            return False
        if t[0]!='OpPhi' or len(t)!=6:return False
        merge=where[value]
        headers=[k for k,lines in blocks.items() if 'OpSelectionMerge '+merge+' None' in lines]
        if len(headers)!=1:return False
        header=headers[0];bs=[x.split() for x in blocks[header] if x.startswith('OpBranchConditional ')]
        if len(bs)!=1 or len(bs[0])!=4 or not condition(bs[0][1]):return False
        branch=bs[0];incoming={t[i+1]:t[i] for i in (2,4)}
        if len(incoming)!=2 or pred[merge]!=set(incoming):return False
        def arm(dest,check):
            rows=[k for k,v in incoming.items() if check(v)]
            if len(rows)!=1:return False
            label=rows[0];value=origin(incoming[label])
            if label==dest:return pred[dest]=={header} and blocks[dest][-1]=='OpBranch '+merge
            nested=choices.get(value)
            return nested is not None and nested['entry']==dest and nested['merge']==label and pred[dest]=={header} and blocks[label][-1]=='OpBranch '+merge
        if not arm(branch[2],yes) or not arm(branch[3],no):return False
        choices[value]={'entry':header,'merge':merge};return True
    def boolconst(value,b):return word(origin(value))[0]==('OpConstantTrue' if b else 'OpConstantFalse')
    def logical(value,op,left,right):
        if binary(value,op,left,right,True):
            choices[origin(value)]={'entry':where[origin(value)],'merge':where[origin(value)]};return True
        return choice(value,left,right,lambda x:boolconst(x,False)) if op=='OpLogicalAnd' else choice(value,left,lambda x:boolconst(x,True),right)
    def true_arm(condition,destination):
        label=where[origin(condition)];b=one(x.split() for x in blocks[label] if x.startswith('OpBranchConditional '))
        assert origin(b[1])==origin(condition) and b[2]==destination and pred[destination]=={label}
    divides=[k for k in local if word(k)[:2]==['OpFDiv',wide]]
    assert len(divides)==2
    rows=[]
    for divide in divides:
        scaled,den=word(divide)[2:];assert word(origin(scaled))[:2]==['OpFMul',wide]
        converted,scale=word(origin(scaled))[2:];assert number(scale,67108864)
        assert word(origin(converted))[:2]==word(origin(den))[:2]==['OpConvertUToF',wide]
        m,d=origin(word(origin(converted))[2]),origin(word(origin(den))[2])
        bits=one(k for k in local if word(k)==['OpBitcast',vec,divide])
        lo=one(k for k in local if word(k)==['OpCompositeExtract',u32,bits,'0'])
        hi=one(k for k in local if word(k)==['OpCompositeExtract',u32,bits,'1'])
        range_cond=one(k for k in local if binary(k,'OpIEqual',lambda v:binary(v,'OpBitwiseAnd',same(hi),num(0xfff00000),True),num(0x41900000),True))
        guess=one(k for k in local if binary(k,'OpBitwiseOr',lambda v:binary(v,'OpShiftLeftLogical',lambda x:binary(x,'OpBitwiseOr',lambda y:binary(y,'OpBitwiseAnd',same(hi),num(0xfffff),True),num(0x100000),True),num(6)),lambda v:binary(v,'OpShiftRightLogical',same(lo),num(26)),True))
        product=one(k for k in local if word(k)[0]=='OpUMulExtended' and [origin(x) for x in word(k)[2:]]==[origin(guess),d])
        assert word(word(product)[1])==['OpTypeStruct',u32,u32]
        def member(index):
            values=[k for k in local if word(k)==['OpCompositeExtract',u32,product,str(index)]]
            assert 1<=len(values)<=2
            for value in values[1:]:aliases[value]=values[0]
            return values[0]
        plo,phi=member(0),member(1)
        nlo=one(k for k in local if binary(k,'OpShiftLeftLogical',same(m),num(26)))
        nhi=one(k for k in local if binary(k,'OpShiftRightLogical',same(m),num(6)))
        below=lambda v:binary(v,'OpULessThan',same(phi),same(nhi))
        equal=lambda v:binary(v,'OpIEqual',same(phi),same(nhi),True)
        lowok=lambda v:binary(v,'OpULessThanEqual',same(plo),same(nlo))
        not_above=one(k for k in local if logical(k,'OpLogicalOr',below,lambda v:logical(v,'OpLogicalAnd',equal,lowok)))
        rlo=one(k for k in local if binary(k,'OpISub',same(nlo),same(plo)))
        borrow_test=one(k for k in local if binary(k,'OpULessThan',same(nlo),same(plo)))
        borrow=lambda v:choice(v,same(borrow_test),num(1),num(0))
        rhi=one(k for k in local if binary(k,'OpISub',lambda v:binary(v,'OpISub',same(nhi),same(phi)),borrow))
        verified=one(k for k in local if logical(k,'OpLogicalAnd',lambda v:binary(v,'OpIEqual',same(rhi),num(0),True),lambda v:binary(v,'OpULessThan',same(rlo),same(d))))
        sticky=lambda v:choice(v,lambda x:binary(x,'OpINotEqual',same(rlo),num(0),True),num(1),num(0))
        jam=one(k for k in local if binary(k,'OpBitwiseOr',same(guess),sticky,True))
        true_arm(range_cond,where[product]);assert where[guess]==where[product]==choices[origin(not_above)]['entry']
        true_arm(not_above,where[rlo]);borrow_value=origin(word(rhi)[3])
        assert choices[borrow_value]['entry']==where[rlo] and choices[borrow_value]['merge']==where[rhi]
        assert choices[origin(verified)]['entry']==where[rhi]
        t=word(jam);sticky_value=origin(t[3] if origin(t[2])==origin(guess) else t[2])
        true_arm(verified,choices[sticky_value]['entry']);assert choices[sticky_value]['merge']==where[jam]
        assert any(binary(k,'OpIEqual',same(jam),num(0),True) for k in local)  # Inlined PB32Round's original zero gate.
        rows.append({'estimate':divide,'normalized_numerator':m,'denominator':d,'bitcast':bits,'range_gate':range_cond,'guess':guess,'exact_product':product,'numerator_words':[nlo,nhi],'remainder_words':[rlo,rhi],'verified':verified,'jam':jam,'closed_guarded_success_and_exact_remainder_footprint':True})
    assert candidate.count('OpUDiv')==stock.count('OpUDiv')==3
    assert candidate.count('OpLoopMerge')==stock.count('OpLoopMerge')
    assert all(word(k)[1]==u32 for k in local if word(k)[0] in ['OpUDiv','OpUMod','OpSRem','OpSMod','OpSDiv'])
    return {'sites':rows,'inlined_camera':True,'integer_division_and_loop_counts_unchanged':True,'normalization_and_complete_common_round_equivalence':'Frozen identical source outside insertion, canonical O1 inlining; not independently reconstructed by this footprint reader.','fallback_scope':'Original complete source fallback exact; emitted integer division/loop counts unchanged. Full named-body fallback comparison applies to the other13 changed modules.'}
