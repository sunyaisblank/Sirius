"""Original function normalization with explicit clone identity and one integer import."""
import re


def functions(text):
    definitions=dict(re.findall(r'^\s*(%\S+) = (.+)$',text,re.M))
    names=dict(re.findall(r'OpName (%\S+) "([^"]*)"',text))
    counts={};identity={}
    for identifier,expression in definitions.items():
        if not expression.startswith('OpFunction '):continue
        name=names[identifier];ordinal=counts.get(name,0);counts[name]=ordinal+1
        identity[identifier]=(name,ordinal)
    return definitions,identity


def normalized_function(text,identifier):
    definitions,identity=functions(text)
    assert identifier in identity
    found=re.search(r'^\s*'+re.escape(identifier)+r' = OpFunction\b.*?OpFunctionEnd',text,re.M|re.S)
    assert found,identifier
    body=found.group(0)+'\n'
    local={value:'local'+str(i) for i,value in enumerate(re.findall(r'^\s*(%\S+) =',body,re.M))}
    def value(identifier):
        if identifier in local:return local[identifier]
        expression=definitions[identifier]
        if expression.startswith('OpFunction '):return ['function',*identity[identifier]]
        if expression.startswith('OpExtInstImport '):
            # The only new external-expression case: exact existing integer FindUMsb import.
            assert expression=='OpExtInstImport "GLSL.std.450"',expression
            return ['OpExtInstImport','GLSL.std.450']
        assert expression.startswith(('OpType','OpConstant','OpSpecConstant')),expression
        return [value(word) if word.startswith('%') else word for word in expression.split()]
    return [[value(word) if word.startswith('%') else word for word in line.split()]
            for line in body.splitlines() if line.strip()]


def compare_helpers(stock,candidate):
    left,left_ids=functions(stock);right,right_ids=functions(candidate)
    left_pb={identity:identifier for identifier,identity in left_ids.items() if identity[0].startswith('PB32')}
    right_pb={identity:identifier for identifier,identity in right_ids.items() if identity[0].startswith('PB32')}
    assert left_pb and set(left_pb)==set(right_pb)
    assert sorted(key for key in left_pb if key[0]=='PB32Divide')==[('PB32Divide',0)]
    compared=[]
    for identity in sorted(left_pb):
        if identity[0]=='PB32Divide':continue
        a,b=left_pb[identity],right_pb[identity]
        assert normalized_function(stock,a)==normalized_function(candidate,b),identity
        compared.append({'name':identity[0],'occurrence':identity[1],'stock_function_id':a,'candidate_function_id':b})
    return {'PB32_function_bodies':len(left_pb),'unique_PB32_names':len({name for name,_ in left_pb}),'unchanged_non_Divide_helper_bodies':len(compared),
            'all_present_helpers':compared,'external_function_reference_identity':'(OpName, definition-order occurrence ordinal)',
            'additional_external_import_supported':'exact OpExtInstImport "GLSL.std.450" only',
            'original_type_constant_specconstant_normalization_rules_retained':True}

