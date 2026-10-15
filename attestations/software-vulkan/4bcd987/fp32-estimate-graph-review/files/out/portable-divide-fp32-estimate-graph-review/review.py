"""Owned pure graph inspection of immutable full Endpoint outputs; no tool/provider invocation."""
from pathlib import Path
import ast,hashlib,json,re,resource,subprocess,sys
sys.dont_write_bytecode=True
from proof import W,R,P,record,content,completed,oldwork,products,archive


def write(path,value):
    path=Path(path);assert not path.exists(),path
    path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')

def static():
    for name,pin in json.loads((W/'source-freeze.json').read_text())['pins'].items():assert record(name)==pin,name

def context():
    c=json.loads((W/'operator-context.json').read_text());assert c['revision']==P['operator_revision']
    for key,name in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:assert c[key]==record(W/name)
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    return c

def definitions():
    stop=Path(P['stop_archive'])/'files'/P['stop_out']
    tree=ast.parse((stop/'screen.py').read_text());names={'function_body','interface_signature'}
    module=ast.Module(body=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names],type_ignores=[])
    assert len(module.body)==2
    scope={'re':re,'json':json};exec(compile(module,'immutable_original_structural_definitions','exec'),scope)
    return scope

def semantic_instructions(text):
    ids=re.findall(r'^\s*(%\S+) =',text,re.M);assert len(ids)==len(set(ids))
    aliases={identifier:'id'+str(i) for i,identifier in enumerate(ids)}
    answer=[]
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith(';'):continue
        tokens=re.findall(r'"(?:[^"\\]|\\.)*"|[^\s]+',line)
        answer.append([aliases[token] if token.startswith('%') else token for token in tokens])
    return answer

def module_checks(proof):
    h=definitions();stock=Path(proof['stock_assembly']['resolved_path']).read_text()
    candidate=(products()/'candidate-endpoint-normal.final.spvasm').read_text()
    assembly=(products()/'candidate-endpoint-normal.spvasm').read_text()
    raw=Path(proof['raw_assembly']['resolved_path']).read_text()
    annotation='\n               OpDecorate %proposal NoContraction'
    location=list(re.finditer(r'^.*OpDecorate .*?$',raw,re.M))[-1].end()
    annotated=raw[:location]+annotation+raw[location:]
    # Reuse the exact already-completed frozen-policy receipt. Do not repeat policy execution.
    entry=re.search(r'OpEntryPoint GLCompute (%\S+)',raw)[1]
    local_size=re.search(r'^.*OpExecutionMode '+re.escape(entry)+r' LocalSize.*$',raw,re.M)[0]
    cap='OpCapability Shader\n               OpCapability RoundingModeRTE'
    rte='\n               OpExecutionMode '+entry+' RoundingModeRTE 32'
    expected=annotated.replace('OpCapability Shader',cap,1).replace(local_size,local_size+rte,1)
    assert assembly==expected and assembly.count(annotation)==1 and assembly.count(rte)==1
    assert assembly.replace(cap,'OpCapability Shader',1).replace(rte,'',1).replace(annotation,'',1)==raw
    # Complete tool-input/readback identity after definition-order ID renaming, not a selected-body subset.
    instructions=semantic_instructions(assembly);assert instructions==semantic_instructions(candidate) and len(instructions)==10339
    assert len(re.findall(r'OpDecorate %proposal NoContraction(?:\s|$)',candidate))==1
    assert sorted(re.findall(r'OpCapability (\S+)',candidate))==['RoundingModeRTE','Shader']
    assert set(re.findall(r'OpTypeInt (\d+) [01]',candidate))=={'32'} and re.findall(r'OpTypeFloat (\d+)',candidate)==['32']
    assert not re.search(r'\b(?:OpConvertFTo[US]|OpConvert[US]ToF|NoUnsignedWrap|NoSignedWrap|FPFastMathMode|FPFastMathDefault|RelaxedPrecision)\b',candidate)
    assert h['interface_signature'](stock)==h['interface_signature'](candidate)
    modes=lambda text:sorted(re.findall(r'OpExecutionMode %\S+ (.*)$',text,re.M))
    assert modes(stock)==modes(candidate)
    def barriers(text):
        defs=dict(re.findall(r'^\s*(%\S+) = (.+)$',text,re.M))
        return [[op,*[defs[word] for word in operands.split()]] for op,operands in re.findall(r'^\s*(OpControlBarrier|OpMemoryBarrier) (.+)$',text,re.M)]
    assert barriers(stock)==barriers(candidate)
    import helper_normalize
    helpers=helper_normalize.compare_helpers(stock,candidate)
    unchanged=helpers['all_present_helpers']
    assert all(any(row['name']==name for row in unchanged) for name in ('PB32LeadingBit','PB32Unpack','PB32Round','PB32Multiply','PB32ProductResidual'))
    wide={'stock':bool(re.search(r'OpName %\S+ "PB32WideProduct"',stock)),'candidate':bool(re.search(r'OpName %\S+ "PB32WideProduct"',candidate))}
    assert wide['stock']==wide['candidate']
    if not wide['stock']:
        assert 'OpUMulExtended' in h['function_body'](stock,'PB32Multiply') and 'OpUMulExtended' in h['function_body'](candidate,'PB32Multiply')
    import graph,controls
    actual=graph.inspect(stock,candidate,h['function_body'])
    tests=controls.run(stock,candidate,h['function_body'],graph.inspect)
    tests['helper_normalization_controls']=helper_normalize.mutation_controls(stock,candidate)
    return {'schema':'sirius-fp32-estimate-full-endpoint-graph-result-v1','passed':True,'actual_Divide':actual,
            'unchanged_normalized_PB32_helper_instances':helpers,'named_PB32WideProduct_presence':wide,
            'inlined_wide_product_semantics':'OpUMulExtended exact UInt2 product; actual helper bodies normalized equal; no standalone function invented.',
            'completed_original_policy_receipt_reused_without_rerun':True,'single_NoContraction_FDiv_scope_retained':True,'assembler_and_readback_complete_instructions_equal':len(instructions),
            'caps':['RoundingModeRTE','Shader'],'all_integer_widths':[32],'all_float_widths':[32],
            'complete_interface_modes_and_ordered_barriers_equal_stock':True,'original_normal_sum_policy_retained':True,
            'products':proof['candidate_products'],'stock_module':proof['stock_module'],'stock_assembly':proof['stock_assembly'],
            'candidate_header_sha256':P['candidate_header_sha256'],'candidate_macro':P['candidate_macro'],
            'native_accuracy_or_RTE_assumption':False,'arithmetic_execution':False,'cost_or_hit_rate_or_adoption_established':False},tests

def pin_context():
    static();completed()
    assert not (W/'operator-context.json').exists() and not (W/'read-graph').exists()
    assert not any((W/name).exists() for name in ('graph-result.json','parser-control-result.json','graph-acceptance.json'))
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()==P['operator_revision']
    assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=R)
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    actions=json.loads((W/'actions.json').read_text());assert list(actions)==['read-graph'] and actions['read-graph']['revision']==P['operator_revision']
    write(W/'operator-context.json',{'revision':P['operator_revision'],'actions':record(W/'actions.json'),'source_freeze':record(W/'source-freeze.json'),'preparation':record(W/'preparation.json')})
    print('Pure graph context pinned; no compiler or provider launched.',flush=True)

def read_graph():
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));_,prior=resource.getrlimit(resource.RLIMIT_FSIZE)
    limit=64*1024**2 if prior==resource.RLIM_INFINITY else min(prior,64*1024**2);resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
    static();context();proof=completed()
    assert not (W/'graph-result.json').exists() and not (W/'parser-control-result.json').exists()
    result,tests=module_checks(proof)
    result['completed_policy_chain']=proof
    static()
    write(W/'parser-control-result.json',tests);write(W/'graph-result.json',result)
    write(W/'read-graph/output-seal.json',{'label':'read-graph','source_freeze':record(W/'source-freeze.json'),'operator_context':record(W/'operator-context.json'),
          'inputs_before':record(W/'read-graph/inputs-before.json'),'outputs':{str(W/name):record(W/name) for name in ('graph-result.json','parser-control-result.json')},
          'arithmetic_execution':False,'cost_or_runtime_or_adoption_admission':False})
    print('Actual full Endpoint graph and source mutation controls accepted; no execution/cost claim.',flush=True)

def review():
    static();context();proof=completed();folder=W/'read-graph'
    owner=json.loads((folder/'owner.json').read_text());action=json.loads((W/'actions.json').read_text())['read-graph']
    assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
    assert owner['source']=={'revision':P['operator_revision'],'status':''} and owner['limits']==P['limits'] and owner['argv']==action['argv']
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['stop_reason'] is None
    assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[] and owner['observed_births'] and all(b['absent'] for b in owner['observed_births'])
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    before=json.loads((folder/'inputs-before.json').read_text())
    for name,pin in before.items():
        path=Path(name);path=R/path if not path.is_absolute() else path
        assert {k:record(path)[k] for k in pin}==pin,name
    seal=json.loads((folder/'output-seal.json').read_text())
    for key,path in [('source_freeze',W/'source-freeze.json'),('operator_context',W/'operator-context.json'),('inputs_before',folder/'inputs-before.json')]:assert seal[key]==record(path)
    expected={str(W/name):record(W/name) for name in ('graph-result.json','parser-control-result.json')};assert seal['outputs']==expected
    assert not seal['arithmetic_execution'] and not seal['cost_or_runtime_or_adoption_admission']
    result=json.loads((W/'graph-result.json').read_text());tests=json.loads((W/'parser-control-result.json').read_text())
    assert result['passed'] and tests['passed'] and tests['positive_count']==3 and tests['negative_count']==38
    assert tests['helper_normalization_controls']['passed'] and tests['helper_normalization_controls']['negative_count']==2
    saved=result['completed_policy_chain'];assert {k:v for k,v in saved.items() if k!='fresh_saved_birth_absence'}=={k:v for k,v in proof.items() if k!='fresh_saved_birth_absence'}
    assert all(row['saved_birth_absent'] for rows in saved['fresh_saved_birth_absence'].values() for row in rows)
    assert not result['arithmetic_execution'] and not result['cost_or_hit_rate_or_adoption_established']
    write(W/'graph-acceptance.json',{'schema':'sirius-fp32-estimate-full-endpoint-graph-acceptance-v1','passed':True,
          'graph_owner':record(folder/'owner.json'),'output_seal':record(folder/'output-seal.json'),'graph_result':record(W/'graph-result.json'),'parser_controls':record(W/'parser-control-result.json'),
          'immutable_policy_chain':proof,'original_failed_passed':False,'original_failed_child_exit':1,
          'structural_graph_accepted':True,'arithmetic_or_cost_or_adoption_admitted':False})
    print('Terminal graph-only acceptance saved; arithmetic and ordinary cost remain separate.',flush=True)

if __name__=='__main__':{'pin-context':pin_context,'read-graph':read_graph,'review':review}[sys.argv[1]]()
