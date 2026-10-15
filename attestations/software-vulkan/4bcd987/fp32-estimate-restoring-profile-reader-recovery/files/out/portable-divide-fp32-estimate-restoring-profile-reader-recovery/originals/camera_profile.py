"""Bounded inlined-camera estimate footprint; frozen-source Round is explicit."""
import collections
import re


def integer_division_signature(text):
    """Exact complete-module original32 opcode/signedness/count signature."""
    definitions = dict(re.findall(r'^\s*(%\S+) = (.+)$', text, re.M))
    types = {key:(int(bits),int(signed)) for key,bits,signed in re.findall(r'(%\S+) = OpTypeInt (\d+) ([01])',text)}
    expected = {'OpUDiv':(32,0),'OpUMod':(32,0),'OpSDiv':(32,1),'OpSRem':(32,1),'OpSMod':(32,1)}
    counts = collections.Counter()
    for definition in definitions.values():
        words = definition.split()
        if words[0] not in expected: continue
        assert len(words) == 4 and types[words[1]] == expected[words[0]], definition
        counts[(words[0],*types[words[1]])] += 1
    assert counts
    return [{'opcode':opcode,'bits':bits,'signed':signed,'count':count}
            for (opcode,bits,signed),count in sorted(counts.items())]


def inspect(stock, candidate):
    definitions = dict(re.findall(r'^\s*(%\S+) = (.+)$', candidate, re.M))
    word = lambda value: definitions[value].split()
    def one(values):
        values = list(values)
        assert len(values) == 1, values
        return values[0]
    u32 = one(k for k,v in definitions.items() if v == 'OpTypeInt 32 0')
    f32 = one(k for k,v in definitions.items() if v == 'OpTypeFloat 32')
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
    assert not re.search(r"\b(?:NoUnsignedWrap|NoSignedWrap|RelaxedPrecision|FPFastMathMode|FPFastMathDefault|OpConvertFTo[SU]|OpConvert[SU]ToF)\b",candidate)
    divides=[k for k in local if word(k)[:2]==['OpFDiv',f32]]
    assert len(divides)==2
    entry=next(iter(blocks))
    dom={label:({entry} if label==entry else set(blocks)) for label in blocks}
    for _ in range(len(blocks)+1):
        changed=False
        for label in blocks:
            if label==entry:continue
            assert pred[label]
            new={label}|set.intersection(*(dom[p] for p in pred[label]))
            if dom[label]!=new:dom[label]=new;changed=True
        if not changed:break
    else:raise AssertionError('camera dominance did not converge')
    rows=[]
    for divide in divides:
        assert re.search(r'OpDecorate '+re.escape(divide)+r' NoContraction(?:\s|$)',candidate)
        numerator_cast,denominator_cast=word(divide)[2:]
        assert word(origin(numerator_cast))[:2]==word(origin(denominator_cast))[:2]==['OpBitcast',f32]
        numerator_bits=origin(word(origin(numerator_cast))[2]);denominator_bits=origin(word(origin(denominator_cast))[2])
        # Identify the denominator by its exact normalized mask operand, not by a name.
        den_parts=word(denominator_bits)[2:]
        den_mask=one(v for v in den_parts if word(origin(v))[0]=='OpBitwiseAnd')
        d=one(v for v in word(origin(den_mask))[2:] if not number(v,0x7fffff))
        selector=one(v for v in word(numerator_bits)[2:] if choice(v,lambda x:word(origin(x))[0]=='OpULessThan',num(0x4d000000),num(0x4c800000)))
        selection=word(origin(selector))
        if selection[0]=='OpSelect':norm=origin(selection[2])
        else:
            start=choices[origin(selector)]['entry'];norm=one(line.split()[1] for line in blocks[start] if line.startswith('OpBranchConditional '))
        assert word(norm)[:2]==['OpULessThan',one(k for k,v in definitions.items() if v=='OpTypeBool')]
        base,normal_d=map(origin,word(norm)[2:]);assert normal_d==origin(d)
        assert binary(numerator_bits,'OpBitwiseOr',same(selector),lambda v:binary(v,'OpBitwiseAnd',same(base),num(0x7fffff),True),True)
        m=one(k for k in local if choice(k,same(norm),lambda v:binary(v,'OpShiftLeftLogical',same(base),num(1)),same(base)))
        bits=one(k for k in local if word(k)==['OpBitcast',u32,divide])
        guard=one(k for k in local if binary(k,'OpIEqual',lambda v:binary(v,'OpBitwiseAnd',same(bits),num(0xff800000),True),num(0x4c800000),True))
        h=one(k for k in local if binary(k,'OpShiftLeftLogical',lambda v:binary(v,'OpBitwiseOr',lambda v:binary(v,'OpBitwiseAnd',same(bits),num(0x7fffff),True),num(0x800000),True),num(3)))
        guess=one(k for k in local if binary(k,'OpISub',same(h),num(16)))
        product=one(k for k in local if word(k)[0]=='OpUMulExtended' and [origin(v) for v in word(k)[2:]]==[origin(guess),origin(d)])
        assert word(word(product)[1])==['OpTypeStruct',u32,u32]
        def member(index):
            values=[k for k in local if word(k)==['OpCompositeExtract',u32,product,str(index)]]
            assert 1<=len(values)<=2
            for value in values[1:]:aliases[value]=values[0]
            return values[0]
        plo,phi=member(0),member(1)
        nlo=one(k for k in local if binary(k,'OpShiftLeftLogical',same(m),num(26)))
        nhi=one(k for k in local if binary(k,'OpShiftRightLogical',same(m),num(6)))
        rlo=one(k for k in local if binary(k,'OpISub',same(nlo),same(plo)))
        borrow_test=one(k for k in local if binary(k,'OpULessThan',same(nlo),same(plo)))
        borrow=one(k for k in local if choice(k,same(borrow_test),num(1),num(0)))
        high_difference=one(k for k in local if binary(k,'OpISub',same(nhi),same(phi)))
        rhi=one(k for k in local if binary(k,'OpISub',same(high_difference),same(borrow)))
        bound=one(k for k in local if binary(k,'OpShiftLeftLogical',same(d),num(5)))
        highzero=lambda v:binary(v,'OpIEqual',same(rhi),num(0),True)
        lowok=lambda v:binary(v,'OpULessThan',same(rlo),same(bound))
        verified=one(k for k in local if logical(k,'OpLogicalAnd',highzero,lowok))
        steps=[];rem=rlo;previous=None
        for shift in (4,3,2,1,0):
            threshold=one(k for k in local if binary(k,'OpShiftLeftLogical',same(d),num(shift))) if shift else d
            test=one(k for k in local if binary(k,'OpUGreaterThanEqual',same(rem),same(threshold)))
            mask=one(k for k in local if choice(k,same(test),num(0xffffffff),num(0)))
            selected=one(k for k in local if binary(k,'OpBitwiseAnd',same(threshold),same(mask),True))
            next_rem=one(k for k in local if binary(k,'OpISub',same(rem),same(selected)))
            bit=one(k for k in local if binary(k,'OpBitwiseAnd',num(1<<shift),same(mask),True))
            if previous is None:
                merged=[k for k in local if binary(k,'OpBitwiseOr',num(0),same(bit),True)]
                assert len(merged)<=1;digit=merged[0] if merged else bit
            else:digit=one(k for k in local if binary(k,'OpBitwiseOr',same(previous),same(bit),True))
            steps.append({'shift':shift,'threshold':threshold,'test':test,'mask':mask,'selected':selected,'input_remainder':rem,'remainder':next_rem,'bit':bit,'digit':digit})
            rem=next_rem;previous=digit
        q=one(k for k in local if binary(k,'OpIAdd',same(guess),same(digit),True))
        sticky=one(k for k in local if choice(k,lambda v:binary(v,'OpINotEqual',same(rem),num(0),True),num(1),num(0)))
        jam=one(k for k in local if binary(k,'OpBitwiseOr',same(q),same(sticky),True))
        # Same guard/control obligations as original inlined profile, now all5 steps.
        true_arm(guard,where[product]);assert where[guess]==where[product]
        assert choices[origin(borrow)]['entry']==where[rlo] and choices[origin(borrow)]['merge']==where[rhi]
        first=choices[origin(steps[0]['mask'])]['entry'];true_arm(verified,first)
        success=where[jam]
        assert where[guard] in dom[success] and where[product] in dom[success] and where[verified] in dom[success]
        for step in steps:
            c=choices[origin(step['mask'])]
            assert c['entry'] in dom[success] and c['merge'] in dom[success]
        assert any(binary(k,'OpIEqual',same(jam),num(0),True) for k in local)  # Original inlined Round zero gate.
        rows.append({'proposal':divide,'normalization_predicate':norm,'base':base,'normalized_m':m,'denominator':d,'guard':guard,'guess':guess,'product':product,'numerator_words':[nlo,nhi],'residual_words':[rlo,rhi],'verified':verified,'steps':steps,'q':q,'sticky':sticky,'jam':jam,'closed_guard_and_restoring_dataflow':True})
    signature=integer_division_signature(candidate)
    assert signature==integer_division_signature(stock)
    assert candidate.count('OpUDiv')==stock.count('OpUDiv')==3
    assert candidate.count('OpLoopMerge')==stock.count('OpLoopMerge')
    return {'inlined_camera':True,'sites':rows,'integer_division_signature':signature,'integer_division_and_loop_counts_unchanged':True,'corrective_UDiv_absent':True,
            'normalization_and_complete_common_round_equivalence':'Exact candidate source differs only by accepted insertion; original O1 inlining. This bounded dataflow reader does not independently reconstruct full inlined normalization/Round or complete fallback.',
            'fallback_scope':'Original source fallback exact, emitted original integer opcode/type/count and loop signatures equal. Named13 modules separately prove complete prefix/fallback/helper bodies.',
            'native_accuracy_or_hitrate_assumed':False,'arithmetic_executed':False}
