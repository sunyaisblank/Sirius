"""Compare actual validated default modules with the guarded private candidates."""
from pathlib import Path
import hashlib
import json
import re
WORK=Path(__file__).resolve().parent

def seal(path):
    raw=path.read_bytes()
    return {'path':str(path),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def module(path):
    text=path.read_text()
    functions={m[1]:m[0] for m in re.finditer(r'^\s*(%\S+) = OpFunction .*?OpFunctionEnd',text,re.M|re.S)}
    names=dict(re.findall(r'OpName (%\S+) "([^"]+)"',text))
    prefix=text[:min(text.index(body) for body in functions.values())]
    definitions=dict(re.findall(r'^\s*(%\S+) = ([^\n]+)',prefix,re.M))
    # Friendly disassembly function IDs distinguish the entry wrapper from its
    # same-named implementation. Require those identities unchanged as well.
    by_name={k.lstrip('%'):(k,v) for k,v in functions.items()}
    return text,functions,names,definitions,by_name

def normalize(body,functions,names,definitions):
    locals={k:f'L{i}' for i,k in enumerate(re.findall(r'^\s*(%\S+) = ',body,re.M))}
    cache={}
    def global_value(key,seen=frozenset()):
        if key in functions:return 'FUNCTION:'+key.lstrip('%')+':'+names[key]
        if key in cache:return cache[key]
        assert key not in seen, ('global definition cycle',key)
        assert key in definitions, ('unknown global',key)
        value=re.sub(r'%[^\s,]+',lambda m:global_value(m[0],seen|{key}),definitions[key])
        # Name storage locations as well as their structural declarations.
        if definitions[key].startswith('OpVariable '):value='GLOBAL:'+names.get(key,key)+':'+value
        cache[key]='('+value+')'
        return cache[key]
    return re.sub(r'%[^\s,]+',lambda m:locals[m[0]] if m[0] in locals else global_value(m[0]),' '.join(body.split()))

result=json.loads((WORK/'compiler-result.json').read_text())
assert result['pass'] and result['device_workloads']==0
before=json.loads((WORK/'compiler-inputs-before.json').read_text())
assert before==json.loads((WORK/'compiler-inputs-after.json').read_text())
commands=json.loads((WORK/'compiler-commands.json').read_text())
assert len(commands)==20
for c in commands:
    assert c['status']=='terminal' and c['exit']==0 and c['leader_reaped'] and not c.get('forced_cleanup')
    for key in ('tool','stdout','stderr'):
        record=c[key];assert seal(Path(record['path']))==record
for record in before:
    assert seal(Path(record['path']))==record
reports=[]
for stage in result['stages']:
    kind=stage['kind'];stem='retained_'+kind.lower()
    base=module(WORK/f'baseline-{stem}.spvasm')
    trial=module(WORK/f'candidate-{stem}.spvasm')
    bt,bf,bn,bd,bb=base;ct,cf,cn,cd,cb=trial
    for name in ('baseline','candidate'):
        for key in ('module','assembly'):
            row=stage['outputs'][name][key]; assert seal(Path(row['path']))==row
    assert set(bb)==set(cb)
    changed=[]
    for name in sorted(bb):
        if normalize(bb[name][1],bf,bn,bd)!=normalize(cb[name][1],cf,cn,cd):changed.append(name)
    assert changed==['RetainedProduct32','RetainedProduct32_0'],changed
    rows=lambda text,op:re.findall(r'^\s*'+op+r' ([^\n]+)',text,re.M)
    assert sorted(rows(ct,'OpCapability'))==sorted(rows(bt,'OpCapability')+['FMAKHR'])
    assert sorted(rows(ct,'OpExtension'))==sorted(rows(bt,'OpExtension')+['"SPV_KHR_fma"'])
    for op in ('OpMemoryModel','OpEntryPoint','OpExecutionMode'):
        assert rows(ct,op)==rows(bt,op),op
    assert not re.search(r'Op(?:TypeFloat|TypeInt) 64|SignedZeroInfNanPreserve|FPFastMathMode|FPFastMathDefault|RelaxedPrecision',ct)
    assert rows(ct,'OpExecutionMode')==['%ComputeMain DenormPreserve 32','%ComputeMain LocalSize 64 1 1','%ComputeMain RoundingModeRTE 32']
    for storage in ('Workgroup','Private'):
        assert re.findall(r'^\s*%\S+ = OpVariable %\S+ '+storage+r'[^\n]*',ct,re.M)==re.findall(r'^\s*%\S+ = OpVariable %\S+ '+storage+r'[^\n]*',bt,re.M)
    abi=lambda text:re.findall(r'^\s*(?:OpDecorate (?:%inputs|%outputs) (?:Binding|DescriptorSet) \d+|OpDecorate %_runtimearr_uint ArrayStride \d+|OpMemberDecorate (?:%StructuredBuffer|%RWStructuredBuffer) 0 Offset \d+)',text,re.M)
    assert abi(ct)==abi(bt) and len(abi(ct))==7
    for op in ('OpFAdd','OpFSub','OpControlBarrier'):
        assert len(re.findall(r'\b'+op+r'\b',ct))==len(re.findall(r'\b'+op+r'\b',bt)),op
    assert len(re.findall(r'\bOpFMul\b',ct))==len(re.findall(r'\bOpFMul\b',bt))+2
    fmas=re.findall(r'(%\S+) = OpFmaKHR %\S+ (%\S+) (%\S+) (%\S+)',ct)
    assert len(fmas)==2 and not re.search(r'\bOpFmaKHR\b',bt)
    decorated=set(re.findall(r'OpDecorate (%\S+) NoContraction',ct))
    assert all(row[0] in decorated for row in fmas)
    assert all(key in decorated for key in re.findall(r'(%\S+) = OpF(?:Add|Sub|Mul) ',ct))
    for name in changed:
        body=cb[name][1]
        assert len(re.findall(r'\bOpFmaKHR\b',body))==1
        call_names=[cn[k] for k in re.findall(r'= OpFunctionCall %\S+ (%\S+)',body)]
        assert call_names==['RPScalarNegate'],call_names
        labels=list(re.finditer(r'^\s*%\S+ = OpLabel\n',body,re.M))
        # Each candidate's final block is the entire original Dekker body. Add
        # the unchanged signature/parameters, then compare every instruction,
        # operand, store and return after structural ID normalization.
        fallback=body[:labels[0].start()]+body[labels[-1].start():]
        assert normalize(fallback,cf,cn,cd)==normalize(bb[name][1],bf,bn,bd),name
    reports.append({'kind':kind,'module':stage['outputs']['candidate']['module'],
                    'functions':len(cf),'only_changed_functions':changed,
                    'complete_original_fallbacks_exact_after_id_normalization':len(changed),
                    'other_functions_exact_after_global_and_local_id_normalization':len(cf)-len(changed),
                    'baseline_execution_modes_entry_and_abi_exact':True,
                    'global_workgroup_private_declarations_exact':True,
                    'required_fma_capability_extension_only':True,'fma_sites':2,
                    'float64_and_relaxed_modes_absent':True,
                    'whole_body_cfg_dominance_and_fallback_requires_independent_review':True})
report={'pass':True,'direct_compiler_inputs_verified':len(before),'stages':reports,
        'scope':'Structural exact module facts; guard operands/CFG reviewed independently. No packed output/device timing/rendering or qualification verdict.'}
(WORK/'compiled-inspection.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
