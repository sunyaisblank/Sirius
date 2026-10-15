"""Bounded actual-FP32 Divide SSA/closed-CFG proof; never invokes a tool."""
import re


def inspect(stock,candidate,function_body):
    module=lambda text:dict(re.findall(r'^\s*(%\S+) = (.+)$',text,re.M))
    definitions=module(candidate);names=dict(re.findall(r'OpName (%\S+) "([^"]*)"',candidate))
    body=function_body(candidate,'PB32Divide');local=module(body)
    def words(value):return definitions[value].split()
    def one(values):
        values=list(values);assert len(values)==1,values;return values[0]
    u32=one(k for k,v in definitions.items() if v=='OpTypeInt 32 0')
    i32=one(k for k,v in definitions.items() if v=='OpTypeInt 32 1')
    boolean=one(k for k,v in definitions.items() if v=='OpTypeBool')
    f32=one(k for k,v in definitions.items() if v=='OpTypeFloat 32')
    assert not re.search(r'\b(?:NoUnsignedWrap|NoSignedWrap)\b',candidate), 'wrapping arithmetic must retain modulo semantics'
    assert len(local)<400 and not re.search(r'\bOp(?:Load|Store|Variable|AccessChain|PtrAccessChain|CopyMemory|CopyObject)\b',body)
    assert len(re.findall(r'= OpUDiv ',body))==2 and body.count('OpLoopMerge')==1
    assert len(re.findall(r'= OpUMulExtended ',body))==1
    assert len(re.findall(r'= OpFDiv ',candidate))==1
    assert re.findall(r'= (OpF(?!unction)\S+) ',body)==['OpFDiv']
    blocks={};location={};current=None
    for raw in body.splitlines():
        line=raw.strip()
        if not line or line=='OpFunctionEnd':continue
        match=re.match(r'(%\S+) = OpLabel$',line)
        if match:current=match[1];blocks[current]=[]
        if current is not None:
            blocks[current].append(line)
            match=re.match(r'(%\S+) =',line)
            if match:location[match[1]]=current
    predecessors={label:set() for label in blocks};successors={label:set() for label in blocks}
    for label,lines in blocks.items():
        row=lines[-1].split();op=row[0]
        assert op in ('OpBranch','OpBranchConditional','OpSwitch','OpReturnValue')
        targets=(row[1:] if op=='OpBranch' else row[2:4] if op=='OpBranchConditional' else [row[2],*row[4::2]] if op=='OpSwitch' else [])
        for target in targets:
            assert target in blocks;predecessors[target].add(label);successors[label].add(target)
    entry=next(iter(blocks));assert not predecessors[entry]
    reachable={entry};pending=[entry]
    while pending:
        for nxt in successors[pending.pop()]:
            if nxt not in reachable:reachable.add(nxt);pending.append(nxt)
    assert reachable==set(blocks)
    dominators={label:({entry} if label==entry else set(blocks)) for label in blocks}
    for _ in range(len(blocks)+1):
        changed=False
        for label in blocks:
            if label==entry:continue
            assert predecessors[label]
            value={label}|set.intersection(*(dominators[p] for p in predecessors[label]))
            if value!=dominators[label]:dominators[label]=value;changed=True
        if not changed:break
    else:raise AssertionError('dominators did not converge')
    def constant(value,number):
        row=words(value);return row[:1]==['OpConstant'] and row[1] in (u32,i32) and int(row[2],0)==number
    number=lambda n:lambda value:constant(value,n)
    same=lambda value:lambda other:value==other
    bool_ops={'OpIEqual','OpINotEqual','OpULessThan','OpULessThanEqual','OpLogicalAnd','OpLogicalOr'}
    def binary(value,opcode,left,right,commute=False):
        row=words(value)
        if len(row)!=4 or row[:2]!=[opcode,boolean if opcode in bool_ops else u32]:return False
        a,b=row[2:];return (left(a) and right(b)) or (commute and left(b) and right(a))
    def isop(value,opcode,a,b):return binary(value,opcode,same(a),same(b))
    choices={}
    def choice(value,condition,yes,no,kind):
        row=words(value)
        if row[0]!='OpPhi' or row[1]!=kind or len(row)!=6:return False
        merge=location[value]
        headers=[label for label,lines in blocks.items() if 'OpSelectionMerge '+merge+' None' in lines]
        if len(headers)!=1:return False
        header=headers[0];branch=blocks[header][-1].split()
        if branch[0]!='OpBranchConditional' or len(branch)!=4 or not condition(branch[1]):return False
        incoming={row[i+1]:row[i] for i in (2,4)}
        if len(incoming)!=2 or predecessors[merge]!=set(incoming):return False
        for dest,predicate in ((branch[2],yes),(branch[3],no)):
            if dest not in incoming or not predicate(incoming[dest]):return False
            if predecessors[dest]!={header} or blocks[dest][-1]!='OpBranch '+merge:return False
        choices[value]={'header':header,'merge':merge,'condition':branch[1],'true':branch[2],'false':branch[3]}
        return True
    def bool_constant(value,truth):return words(value)==['OpConstantTrue' if truth else 'OpConstantFalse',boolean]
    def logical_and(value,left,right):
        if binary(value,'OpLogicalAnd',left,right,True):return True
        return choice(value,left,right,lambda v:bool_constant(v,False),boolean)
    def normalized(text,fragment,roles=None):
        defs=module(text);fnames=dict(re.findall(r'OpName (%\S+) "([^"]*)"',text))
        ids=re.findall(r'^\s*(%\S+) =',fragment,re.M);mapped={value:'local'+str(i) for i,value in enumerate(ids)};roles=roles or {}
        def value(identifier):
            if identifier in mapped:return mapped[identifier]
            if identifier in roles:return roles[identifier]
            row=defs[identifier].split()
            if row[0]=='OpFunction':return ['function',fnames[identifier]]
            assert row[0].startswith(('OpType','OpConstant')) or row[0]=='OpISub',row
            return [value(w) if w.startswith('%') else w for w in row]
        return [[value(w) if w.startswith('%') else w for w in line.split()] for line in fragment.splitlines() if line.strip()]
    def fallback(text):
        defs=module(text);fnames=dict(re.findall(r'OpName (%\S+) "([^"]*)"',text));fragment=function_body(text,'PB32Divide');inside=module(fragment)
        row=lambda k:defs[k].split()
        divisions=[k for k in inside if row(k)[0]=='OpUDiv' and row(row(k)[2])[0]=='OpShiftLeftLogical' and row(row(row(k)[2])[2])[0]=='OpPhi']
        division=one(divisions);num,d=row(division)[2:];rem=row(num)[2]
        loop=one(re.finditer(r'OpLoopMerge (%\S+) (%\S+) (\S+)',fragment))
        start=list(re.finditer(r'^\s*(%\S+) = OpLabel$',fragment[:loop.start()],re.M))[-1].start();suffix=fragment[start:]
        mapped=set(re.findall(r'^\s*(%\S+) =',suffix,re.M));pairs=row(rem)[2:]
        initial_pairs=[(pairs[i],pairs[i+1]) for i in range(0,len(pairs),2) if pairs[i+1] not in mapped];initial,preheader=one(initial_pairs)
        assert row(initial)[:1]==['OpISub'] and row(initial)[3]==d
        m=row(initial)[2];assert row(m)[0]=='OpPhi'
        calls=[row(k) for k in inside if row(k)[0]=='OpFunctionCall' and fnames.get(row(k)[2])=='PB32Round'];tail=one(c for c in calls if c[-1] in mapped)
        roles={m:'normalized_m',d:'divisor',tail[3]:'sign',tail[4]:'exponent',preheader:'single_preheader'}
        external_labels={v for v in re.findall(r'%[^\s]+',suffix) if v in defs and row(v)[0]=='OpLabel' and v not in mapped};assert external_labels=={preheader}
        line=re.search(r'^\s*'+re.escape(m)+r' = .*$',fragment,re.M)
        prefix=fragment[:line.end()]+'\n'
        return {'normalized':normalized(text,suffix,roles),'m':m,'d':d,'sign':tail[3],'exponent':tail[4],'preheader':preheader,'initial':initial,'tail':tail,'prefix':prefix,'division':division,'loop_label':one(k for k in mapped if defs[k]=='OpLabel' and suffix.lstrip().startswith(k+' = OpLabel'))}
    before=fallback(stock);after=fallback(candidate)
    assert before['normalized']==after['normalized']
    assert normalized(stock,before['prefix'])==normalized(candidate,after['prefix'])
    prefix_ids=lambda info:re.findall(r'^\s*(%\S+) =',info['prefix'],re.M)
    for role in ('m','d','sign','exponent'):
        assert prefix_ids(before).index(before[role])==prefix_ids(after).index(after[role]),role
    m,d,sign,exponent=(after[k] for k in ('m','d','sign','exponent'))
    assert words(m)[:2]==['OpPhi',u32]
    mrow=words(m);base_candidates=[mrow[i] for i in (2,4) if words(mrow[i])[:1]==['OpCompositeExtract']]
    base=one(base_candidates);assert words(base)[1]==u32 and words(base)[3:]==['0']
    norm_pred=one(k for k in local if isop(k,'OpULessThan',base,d))
    assert choice(m,same(norm_pred),lambda v:binary(v,'OpShiftLeftLogical',same(base),number(1)),same(base),u32)
    selector=one(k for k in local if choice(k,same(norm_pred),number(0x4d000000),number(0x4c800000),u32))
    numerator_bits=one(k for k in local if binary(k,'OpBitwiseOr',same(selector),lambda v:binary(v,'OpBitwiseAnd',same(base),number(0x7fffff),True),True))
    denominator_bits=one(k for k in local if binary(k,'OpBitwiseOr',number(0x3f800000),lambda v:binary(v,'OpBitwiseAnd',same(d),number(0x7fffff),True),True))
    cast=lambda value,kind,source:words(value)==['OpBitcast',kind,source]
    native_div=one(k for k in local if words(k)[:2]==['OpFDiv',f32])
    assert cast(words(native_div)[2],f32,numerator_bits) and cast(words(native_div)[3],f32,denominator_bits)
    assert re.search(r'OpDecorate '+re.escape(native_div)+r' NoContraction(?:\s|$)',candidate)
    bits=one(k for k in local if cast(k,u32,native_div))
    guard=one(k for k in local if binary(k,'OpIEqual',lambda v:binary(v,'OpBitwiseAnd',same(bits),number(0xff800000),True),number(0x4c800000),True))
    h=one(k for k in local if binary(k,'OpShiftLeftLogical',lambda v:binary(v,'OpBitwiseOr',lambda w:binary(w,'OpBitwiseAnd',same(bits),number(0x7fffff),True),number(0x800000),True),number(3)))
    guess=one(k for k in local if binary(k,'OpISub',same(h),number(16)))
    product=one(k for k in local if words(k)[0]=='OpUMulExtended');assert words(product)[2:]==[guess,d]
    assert words(words(product)[1])==['OpTypeStruct',u32,u32]
    member=lambda i:one(k for k in local if words(k)==['OpCompositeExtract',u32,product,str(i)])
    plo,phi=member(0),member(1)
    nlo=one(k for k in local if binary(k,'OpShiftLeftLogical',same(m),number(26)))
    nhi=one(k for k in local if binary(k,'OpShiftRightLogical',same(m),number(6)))
    rlo=one(k for k in local if isop(k,'OpISub',nlo,plo))
    borrow_test=one(k for k in local if isop(k,'OpULessThan',nlo,plo))
    borrow=one(k for k in local if choice(k,same(borrow_test),number(1),number(0),u32))
    high_difference=one(k for k in local if isop(k,'OpISub',nhi,phi))
    rhi=one(k for k in local if isop(k,'OpISub',high_difference,borrow))
    high_zero=one(k for k in local if binary(k,'OpIEqual',same(rhi),number(0),True))
    bound=one(k for k in local if binary(k,'OpShiftLeftLogical',same(d),number(5)))
    low_small=one(k for k in local if isop(k,'OpULessThan',rlo,bound))
    verified=one(k for k in local if logical_and(k,same(high_zero),same(low_small)))
    digit=one(k for k in local if isop(k,'OpUDiv',rlo,d))
    low_product=one(k for k in local if binary(k,'OpIMul',same(digit),same(d),True))
    rem=one(k for k in local if isop(k,'OpISub',rlo,low_product))
    quotient=one(k for k in local if binary(k,'OpIAdd',same(guess),same(digit),True))
    rem_nonzero=one(k for k in local if binary(k,'OpINotEqual',same(rem),number(0),True))
    sticky=one(k for k in local if choice(k,same(rem_nonzero),number(1),number(0),u32))
    jam=one(k for k in local if binary(k,'OpBitwiseOr',same(quotient),same(sticky),True))
    rounded=one(k for k in local if words(k)==['OpFunctionCall',u32,one(v for v in definitions if names.get(v)=='PB32Round'),sign,exponent,jam])
    def branch(condition):
        row=blocks[location[condition]][-1].split();assert row[:2]==['OpBranchConditional',condition] and len(row)==4;return row
    range_branch=branch(guard);assert range_branch[2]==location[product] and predecessors[location[product]]=={location[guard]}
    assert range_branch[3]==after['preheader']
    assert location[h]==location[guess]==location[product]==location[rlo]==location[borrow_test]==location[high_difference]
    assert choices[borrow]['header']==location[product] and choices[borrow]['merge']==location[rhi]
    assert choices[verified]['header']==location[rhi] and choices[verified]['condition']==high_zero
    assert choices[verified]['true']==location[low_small]==location[bound]
    verified_branch=branch(verified);assert verified_branch[2]==location[digit] and predecessors[location[digit]]=={location[verified]}
    rejected=verified_branch[3];assert predecessors[rejected]=={location[verified]} and blocks[rejected]==[rejected+' = OpLabel','OpBranch '+after['preheader']]
    assert predecessors[after['preheader']]=={location[guard],rejected}
    assert location[digit]==location[low_product]==location[rem]==location[quotient]==choices[sticky]['header']
    assert choices[sticky]['merge']==location[jam]==location[rounded]
    assert blocks[location[rounded]][-1]=='OpReturnValue '+rounded
    assert location[guard] in dominators[location[rounded]] and location[product] in dominators[location[rounded]] and location[verified] in dominators[location[rounded]]
    assert choices[selector]['header']==location[m] and choices[selector]['merge']==location[native_div]==location[guard]
    assert choices[m]['merge']==choices[selector]['header'] and location[native_div] in dominators[location[product]]
    # Exact new-region coverage: matched definitions and terminators only, in their actual order.
    # The old prefix through m and the complete fallback loop/tail are normalized independently.
    rowline=lambda identifier:identifier+' = '+definitions[identifier]
    def exact(label,ids,tail):
        assert blocks[label]==[label+' = OpLabel',*[rowline(identifier) for identifier in ids],*tail],label
    def empty_choice(value):
        c=choices[value]
        for label in (c['true'],c['false']):exact(label,[],['OpBranch '+c['merge']])
    def select_tail(value):
        c=choices[value];return ['OpSelectionMerge '+c['merge']+' None','OpBranchConditional '+c['condition']+' '+c['true']+' '+c['false']]
    norm_block=blocks[location[m]];m_position=norm_block.index(rowline(m))
    assert norm_block[m_position+1:]==select_tail(selector)
    empty_choice(selector)
    numerator_mask=words(numerator_bits)[3];denominator_mask=words(denominator_bits)[3]
    numerator_cast,denominator_cast=words(native_div)[2:]
    range_mask=words(guard)[2]
    exact(location[selector],[selector,numerator_mask,numerator_bits,denominator_mask,denominator_bits,numerator_cast,denominator_cast,native_div,bits,range_mask,guard],['OpSelectionMerge '+after['preheader']+' None','OpBranchConditional '+guard+' '+location[product]+' '+after['preheader']])
    h_or=words(h)[2];h_mask=words(h_or)[2]
    exact(location[product],[h_mask,h_or,h,guess,product,nlo,nhi,plo,rlo,phi,high_difference,borrow_test],select_tail(borrow))
    empty_choice(borrow)
    exact(location[rhi],[borrow,rhi,high_zero],select_tail(verified))
    exact(choices[verified]['true'],[bound,low_small],['OpBranch '+choices[verified]['merge']])
    exact(choices[verified]['false'],[],['OpBranch '+choices[verified]['merge']])
    exact(location[verified],[verified],['OpSelectionMerge '+rejected+' None','OpBranchConditional '+verified+' '+location[digit]+' '+rejected])
    exact(location[digit],[digit,low_product,rem,quotient,rem_nonzero],select_tail(sticky))
    empty_choice(sticky)
    exact(location[jam],[sticky,jam,rounded],['OpReturnValue '+rounded])
    exact(rejected,[],['OpBranch '+after['preheader']])
    exact(after['preheader'],[after['initial']],['OpBranch '+after['loop_label']])
    new_labels={choices[selector]['true'],choices[selector]['false'],location[selector],location[product],choices[borrow]['true'],choices[borrow]['false'],location[rhi],choices[verified]['true'],choices[verified]['false'],location[verified],location[digit],choices[sticky]['true'],choices[sticky]['false'],location[jam],rejected,after['preheader']}
    first_new=body.index(rowline(m))+len(rowline(m));end_new=body.index(after['loop_label']+' = OpLabel')
    region=body[first_new:end_new]
    assert set(re.findall(r'(%\S+) = OpLabel',region))==new_labels
    assert not re.search(r'\b(?:NoUnsignedWrap|NoSignedWrap)\b',candidate)
    return {'native_accuracy_or_RTE_assumed':False,'native_division':native_div,'normalization_predicate':norm_pred,'base':base,'normalized_m':m,'divisor':d,'sign':sign,'exponent':exponent,'binade_selector':selector,'numerator_bits':numerator_bits,'denominator_bits':denominator_bits,'bits':bits,'range_guard':guard,'h':h,'guess':guess,'product':product,'product_lo':plo,'product_hi':phi,'numerator_lo':nlo,'numerator_hi':nhi,'residual_lo':rlo,'borrow_test':borrow_test,'borrow':borrow,'high_difference':high_difference,'residual_hi':rhi,'high_zero':high_zero,'bound':bound,'low_small':low_small,'verified':verified,'digit':digit,'low_product':low_product,'remainder':rem,'quotient':quotient,'remainder_nonzero':rem_nonzero,'sticky':sticky,'jam':jam,'rounded':rounded,'fallback_division':after['division'],'fallback_preheader':after['preheader'],'closed_choices':choices,'guard_and_success_dominators':sorted(dominators[location[rounded]]),'original_complete_prefix_normalized_equal':True,'original_complete_fallback_suffix_normalized_equal':True,'direct_SSA_no_load_store_alias_admission':True,'single_existing_UMulExtended_semantics':True,'exact_new_region_instruction_and_block_coverage':True,'NoUnsignedWrap_and_NoSignedWrap_absent':True,'original_common_round_sign_exponent_join':True,'arithmetic_executed':False,'cost_or_hit_rate_established':False}
