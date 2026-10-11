"""Inspect compact canonical table reads and retained original generic fallback."""
from pathlib import Path
import hashlib, json, re
import inspect_module as helpers
WORK = Path(__file__).resolve().parent

def main():
    result = json.loads((WORK/'compile-result.json').read_text())
    before = json.loads((WORK/'compile-before.json').read_text())
    commands = json.loads((WORK/'commands.json').read_text())
    assert result['pass'] and result['input_checks']['changed'] == []
    assert helpers.seal(WORK/'literal-transport.spv') == result['module']
    assert helpers.seal(WORK/'literal-transport.spvasm') == result['assembly']
    refs = [before[n] for n in ('source','emitter','runner','program_metadata','emission_metadata','schedule','header','builder','baseline','baseline_copy')] + before['stage_inputs'] + list(before['tools'].values())
    assert len(commands) == 7
    for command in commands:
        cleanup = command['cleanup']
        assert command['disposition'] == 'completed' and command['returncode'] == 0 and cleanup['known_leader_reaped'] and not cleanup['forced_cleanup'] and not cleanup['errors'] and cleanup['terminal_group_scan']['members'] == []
        refs += [command['tool'],command['stdout_seal'],command['stderr_seal']]
    for item in refs:
        data=Path(item['path']).read_bytes()
        assert len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['sha256']
    text,functions,constants,names = helpers.read(WORK/'literal-transport.spvasm')
    assert not any(re.fullmatch(r'Literal(?:General|Schwarzschild)(?:Layer\d{3}|Program)', names.get(key,'')) for key in functions)
    _,recognize = helpers.body_by_name(functions,names,'RecognizeLiteralTransport')
    _,main_body = helpers.body_by_name(functions,names,'ComputeMain')
    main_calls = [names.get(key,key) for key in re.findall(r'= OpFunctionCall %\S+ (%\S+)',main_body)]
    assert main_calls.count('RecognizeLiteralTransport') == 1 and main_calls.count('EvaluateProgramGeneric') == 1 and main_calls.count('EvaluateProgramCanonical') == 1
    assert len(re.findall(r'\bOpControlBarrier\b',recognize)) == 2 and len(re.findall(r'\bOpAtomicAnd\b',recognize)) == 2
    assert 'OpVariable' not in recognize
    programs = json.loads((WORK/'programs.json').read_text())
    composites = {key:words.split() for key,words in re.findall(r'^\s*(%\S+) = OpConstantComposite %\S+ ([^\n]+)',text,re.M)}
    tables={}
    for family,program in programs.items():
        expected=[program['instructions'],program['registers'],len(program['outputs']),*program['outputs'],*program['operations'],*program['layer_offsets']]
        matches=[key for key,words in composites.items() if len(words)==len(expected) and all(word in constants for word in words) and [constants[word] for word in words]==expected]
        assert len(matches)==1
        table='%Literal'+family+'Table'
        assert len(re.findall(r'OpStore '+re.escape(table)+r' ',text))==1
        assert len(re.findall(r'Op(?:InBounds)?AccessChain %\S+ '+re.escape(table)+r' ',text))==2
        assert re.search(r'Op(?:InBounds)?AccessChain %\S+ '+re.escape(table)+r' ',recognize)
        tables[family]={'words':len(expected),'exact_constant_composite':matches[0],'private_table':table}
    _,getter=helpers.body_by_name(functions,names,'CanonicalTransportWord')
    _,node=helpers.body_by_name(functions,names,'EvaluateOneCanonical')
    _,canonical=helpers.body_by_name(functions,names,'EvaluateProgramCanonical')
    _,original_node=helpers.body_by_name(functions,names,'EvaluateOne')
    _,generic=helpers.body_by_name(functions,names,'EvaluateProgramGeneric')
    calls=lambda body:[names.get(key,key) for key in re.findall(r'= OpFunctionCall %\S+ (%\S+)',body)]
    assert len(functions)==96 and 'OpVariable' not in getter and 'OpStore' not in getter
    assert len(re.findall(r'Op(?:InBounds)?AccessChain %\S+ %Literal(?:General|Schwarzschild)Table ',getter))==2
    assert calls(node).count('CanonicalTransportWord')==5 and calls(canonical).count('CanonicalTransportWord')==6
    assert [c for c in calls(node) if c!='CanonicalTransportWord']==calls(original_node)
    assert [c.replace('EvaluateOneCanonical','EvaluateOne') for c in calls(canonical) if c!='CanonicalTransportWord']==calls(generic)
    assert len(re.findall(r'\bOpControlBarrier\b',canonical))==len(re.findall(r'\bOpControlBarrier\b',generic))==4
    baseline_text,_,_,_=helpers.read(WORK/'baseline-transport.spvasm')
    for op in ('OpCapability','OpExtension','OpMemoryModel','OpExecutionMode'):
        assert re.findall(r'^\s*'+op+r' [^\n]+',text,re.M)==re.findall(r'^\s*'+op+r' [^\n]+',baseline_text,re.M),op
    entry=lambda module:re.findall(r'^\s*OpEntryPoint ([^\n]+)',module,re.M)[0].split()
    current_entry,base_entry=entry(text),entry(baseline_text)
    assert current_entry[:3]==base_entry[:3] and set(current_entry[3:])-set(base_entry[3:])=={'%literalTransportMatch','%LiteralSchwarzschildTable','%LiteralGeneralTable'} and not set(base_entry[3:])-set(current_entry[3:])
    abi=lambda module:re.findall(r'^\s*(?:OpDecorate (?:%inputs|%outputs) (?:Binding|DescriptorSet) \d+|OpDecorate %_runtimearr_uint ArrayStride \d+|OpMemberDecorate (?:%StructuredBuffer|%RWStructuredBuffer) 0 Offset \d+)',module,re.M)
    assert abi(text)==abi(baseline_text)
    for op in ('OpFAdd','OpFSub','NoContraction'):
        assert len(re.findall(r'\b'+op+r'\b',text))==len(re.findall(r'\b'+op+r'\b',baseline_text)),op
    assert json.loads((WORK/'source-frame-check.json').read_text())['whole_source_inverse_byte_exact']
    storage=helpers.storage_facts(text,constants)
    assert storage['declared_workgroup_member_bytes']==9688
    assert sum(v['logical_member_bytes'] for v in storage['variables'] if v['id'] in ('%LiteralGeneralTable','%LiteralSchwarzschildTable'))==63868
    report={'pass':True,'scope':'Compact canonical table-read compiler screen; original generic fallback retained. No device/numerical/benefit result. Main CFG, loop/composite-store joins require independent review.','module':result['module'],'assembly':result['assembly'],'seals_checked':len(refs),'function_count':len(functions),'literal_functions_present':0,'recognizer_entry_calls':1,'generic_stage_calls':1,'canonical_stage_calls':1,'canonical_node_table_reads':5,'canonical_program_table_reads':6,'canonical_barriers':4,'canonical_primitive_call_sequence_exact_original':True,'recognizer_barriers':2,'recognizer_atomic_mismatch_sites':2,'complete_tables':tables,'storage':storage,'main_calls':main_calls}
    (WORK/'compiled-inspection.json').write_text(json.dumps(report,indent=2)+'\n')
    for name,body in [('RecognizeLiteralTransport',recognize),('ComputeMain',main_body),('CanonicalTransportWord',getter),('EvaluateOneCanonical',node),('EvaluateProgramCanonical',canonical),('EvaluateProgramGeneric',generic)]:
        (WORK/(name+'.spvasm')).write_text(body+'\n')
    print(json.dumps({key:report[key] for key in ('pass','function_count','literal_functions_present','recognizer_entry_calls','generic_stage_calls','seals_checked')}))

if __name__=='__main__':
    main()
