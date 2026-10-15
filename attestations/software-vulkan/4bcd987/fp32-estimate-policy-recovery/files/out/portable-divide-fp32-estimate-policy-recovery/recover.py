"""Finish one immutable FP32 compiler-policy stop; no compilation or provider."""
from pathlib import Path
import ast,hashlib,importlib.util,json,os,re,resource,stat,subprocess,sys
sys.dont_write_bytecode=True
W=Path(__file__).resolve().parent
R=W.parents[1]
B=R/'bin/linux-gcc/portable-divide-fp32-estimate-policy-recovery'
P=json.loads((W/'packet.json').read_text())

def record(path):
    path=Path(path)
    with path.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':h,'mode':stat.S_IMODE(path.stat().st_mode),
            'resolved_path':str(path.resolve(strict=True)),'symlink_target':os.readlink(path) if path.is_symlink() else None}

def content(path):
    v=record(path);return {k:v[k] for k in ('bytes','sha256')}

def write(path,value):
    path=Path(path);assert not path.exists(),path
    path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')

def static():
    for name,pin in json.loads((W/'source-freeze.json').read_text())['pins'].items():assert record(name)==pin,name

def completed():
    assert re.fullmatch('[0-9a-f]{40}',P['stop_commit']), 'immutable stop commit pending'
    assert re.fullmatch('[0-9a-f]{64}',P['stop_manifest_sha256']), 'immutable stop manifest pending'
    archive=Path(P['stop_archive']);manifest=archive/'manifest.json'
    assert content(manifest)['sha256']==P['stop_manifest_sha256']
    for name,pin in json.loads(manifest.read_text())['files'].items():
        path=archive/name;assert not path.is_symlink() and content(path)==pin,name
    old=archive/'files'/P['original_out'];products=archive/'files'/P['original_bin']/'compiler/modules'
    def joined(name,pin):
        original=Path(name);original=R/original if not original.is_absolute() else original
        owned=any(original.is_relative_to(R/P[key]) for key in ('original_out','original_bin'))
        actual=archive/'files'/original.relative_to(R) if owned else original
        assert content(actual)=={k:pin[k] for k in ('bytes','sha256')},name
        assert pin['resolved_path']==str(original if owned else actual.resolve(strict=True)),name
        if 'mode' in pin:assert record(actual)['mode']==pin['mode']
        if 'symlink_target' in pin:assert record(actual)['symlink_target']==pin['symlink_target']
        return actual
    frozen=json.loads((old/'source-freeze.json').read_text());assert len(frozen['private_inputs'])+len(frozen['external_inputs'])==96
    for name,pin in frozen['private_inputs'].items():joined(R/P['original_out']/name,pin)
    for name,pin in frozen['external_inputs'].items():joined(name,pin)
    context=json.loads((old/'operator-context.json').read_text())
    assert context['revision']==P['operator_revision'] and context['status']==''
    for key,name in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:joined(R/P['original_out']/name,context[key])
    folder=old/'compile-endpoint-normal';owner=json.loads((folder/'owner.json').read_text())
    assert content(folder/'owner.json')==P['failed_owner']
    assert owner['status']=='terminal' and owner['passed'] is False and owner['child_exit']==1
    assert owner['source']=={'revision':P['operator_revision'],'status':''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['stop_reason'] is None
    assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert len(owner['observed_births'])==2 and all(b['absent'] for b in owner['observed_births'])
    assert owner['limits']==P['limits']
    assert owner['argv']==['/usr/bin/python3','-B',str(R/P['original_out']/'screen.py'),'compile']
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    before=json.loads((folder/'inputs-before.json').read_text());assert len(before)==630
    for name,pin in before.items():joined(name,pin)
    stderr=(folder/'stderr.log').read_text()
    assert stderr.endswith('ValueError: diagnostic FP32 proposal lost its single scoped division/order\n')
    assert 'diagnostic_divide_fp32_normal_sum_controls' in stderr
    assert not (old/'compiler-products.json').exists() and not (folder/'output-seal.json').exists()
    assert not (products/'candidate-endpoint-normal.spv').exists()
    assert content(products/'candidate-endpoint-normal.compiler.spv')==P['raw_module']
    assert content(products/'candidate-endpoint-normal.spvasm')==P['raw_assembly']
    builder=old/'candidate/scripts/build-retained-kernels.py'
    assert content(builder)['sha256']==P['frozen_builder_sha256']
    assert content(old/'candidate/src/sirius/kernels/portable_binary32.h')['sha256']==P['candidate_header_sha256']
    independent=old/'independent-compiler-policy-stop-review.json'
    assert content(independent)['sha256']==P['stop_independent_review_sha256'] and json.loads(independent.read_text())['passed']
    reuse=json.loads((old/'stock-reuse.json').read_text());history=Path(reuse['historical_work'])
    historical=json.loads((history/'compile-endpoint-normal/inputs-before.json').read_text())
    history_owner=json.loads((history/'compile-endpoint-normal/owner.json').read_text())
    assert history_owner['status']=='terminal' and history_owner['passed'] and history_owner['child_exit']==0
    assert history_owner['stop_reason'] is None and history_owner['cleanup_errors']==[] and history_owner['remaining_owned_processes']==[]
    assert history_owner['source_unchanged'] and history_owner['whole_inputs_unchanged'] and all(b['absent'] for b in history_owner['observed_births'])
    assert (history/'compile-endpoint-normal/inputs-before.json').read_bytes()==(history/'compile-endpoint-normal/inputs-after.json').read_bytes()
    for row in reuse['source_joins']:
        expected={k:historical[row['original_input_key']][k] for k in ('bytes','sha256')}
        for key in ('canonical','historical','private'):
            original=Path(row[key])
            actual=old/original.relative_to(R/P['original_out']) if key=='private' else original
            assert content(actual)==expected,(key,original)
    for row in reuse['tool_joins']:
        pin=historical[row['original_input_key']];assert content(row['path'])=={k:pin[k] for k in ('bytes','sha256')}
        assert str(Path(row['path']).resolve(strict=True))==pin['resolved_path']
    args=json.loads((old/'generated/endpoint-arguments.json').read_text())
    assert (args['registers'],args['terms'],args['layers'],args['prefix'])==(612,5,111,250)
    assert args['endpoint_specialized']==[11240,9903,84,202] and args['tools']==P['tools']
    assert (612*5+7+16*5)*4+8==12596<16384
    assert content(reuse['stock_module'])=={'bytes':151400,'sha256':'981bc85d13e168c6adf72c99adfe169157593dda1bee46e78ff867fc6b199ede'}
    return {'manifest':record(manifest),'old_failed_owner':record(folder/'owner.json'),
            'old_failed_passed':False,'old_failed_child_exit':1,'raw_module':record(products/'candidate-endpoint-normal.compiler.spv'),
            'raw_assembly':record(products/'candidate-endpoint-normal.spvasm'),'independent_stop_review':record(independent),
            'frozen_builder':record(builder),'source_freeze':record(old/'source-freeze.json'),
            'old630_consumed_inputs_current_or_explicit_archive_equal':True,
            'stock_module':record(reuse['stock_module']),'stock_assembly':record(reuse['stock_assembly'])}

def helpers(old):
    # Definitions only: the frozen producer's main/footer is never run.
    tree=ast.parse((old/'screen.py').read_text());names={'function_body','interface_signature'}
    module=ast.Module(body=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names],type_ignores=[])
    scope={'re':re,'json':json};exec(compile(module,'frozen_interface_definitions','exec'),scope)
    return scope

def pin_context():
    static();completed()
    assert not B.exists() and not (W/'operator-context.json').exists() and not (W/'finish-policy').exists()
    assert not any((W/n).exists() for n in ('policy-result.json','policy-acceptance.json'))
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()==P['operator_revision']
    assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=R)
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    a=json.loads((W/'actions.json').read_text());assert list(a)==['finish-policy'] and a['finish-policy']['revision']==P['operator_revision']
    write(W/'operator-context.json',{'revision':P['operator_revision'],'actions':record(W/'actions.json'),'source_freeze':record(W/'source-freeze.json'),'preparation':record(W/'preparation.json')})
    print('Recovery context pinned; no compiler or provider launched.',flush=True)

def finish():
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));_,prior=resource.getrlimit(resource.RLIMIT_FSIZE)
    limit=64*1024**2 if prior==resource.RLIM_INFINITY else min(prior,64*1024**2);resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
    static();proof=completed();assert not B.exists() and not (W/'policy-result.json').exists()
    context=json.loads((W/'operator-context.json').read_text())
    assert context['revision']==P['operator_revision'] and context['source_freeze']==record(W/'source-freeze.json')
    assert context['actions']==record(W/'actions.json') and context['preparation']==record(W/'preparation.json')
    old=Path(P['stop_archive'])/'files'/P['original_out'];h=helpers(old)
    raw=Path(proof['raw_assembly']['resolved_path']).read_text();body=h['function_body'](raw,'PB32Divide')
    floating=re.findall(r'(%\S+) = OpTypeFloat (\d+)',raw);assert floating==[('%float','32')]
    utypes=re.findall(r'OpTypeInt (\d+) [01]',raw);assert utypes and set(utypes)=={'32'}
    assert re.findall(r'OpCapability (\S+)',raw)==['Shader']
    divisions=re.findall(r'^\s*(%\S+) = OpFDiv %float .*?$',raw,re.M)
    assert divisions==['%proposal'] and re.findall(r'= (OpF(?!unction)\S+) ',body)==['OpFDiv']
    assert re.search(r'^\s*%proposal = OpFDiv %float ',body,re.M)
    assert not re.search(r'^\s*OpDecorate %proposal\b',raw,re.M)
    annotation='\n               OpDecorate %proposal NoContraction'
    location=list(re.finditer(r'^.*OpDecorate .*?$',raw,re.M))[-1].end()
    annotated=raw[:location]+annotation+raw[location:]
    assert annotated.replace(annotation,'',1)==raw
    spec=importlib.util.spec_from_file_location('unchanged_frozen_policy',old/'candidate/scripts/build-retained-kernels.py')
    policy=importlib.util.module_from_spec(spec);spec.loader.exec_module(policy)
    transformed=policy.diagnostic_divide_fp32_normal_sum_controls(annotated)
    entry=re.search(r'OpEntryPoint GLCompute (%\S+)',raw)[1]
    local_size=re.search(r'^.*OpExecutionMode '+re.escape(entry)+r' LocalSize.*$',raw,re.M)[0]
    cap='OpCapability Shader\n               OpCapability RoundingModeRTE'
    rte='\n               OpExecutionMode '+entry+' RoundingModeRTE 32'
    expected=annotated.replace('OpCapability Shader',cap,1).replace(local_size,local_size+rte,1)
    assert transformed==expected
    assert transformed.count(annotation)==1 and transformed.count(rte)==1
    assert transformed.replace(cap,'OpCapability Shader',1).replace(rte,'',1).replace(annotation,'',1)==raw
    assert sorted(re.findall(r'OpCapability (\S+)',transformed))==['RoundingModeRTE','Shader']
    assert 'OpExecutionMode '+entry+' LocalSize 16 1 1' in transformed and 'OpControlBarrier' in transformed
    # Remaining original compile_shader tail: complete full assembly, Vulkan1.2 validation, readback.
    B.mkdir();assembly=B/'candidate-endpoint-normal.spvasm';module=B/'candidate-endpoint-normal.spv';readback=B/'candidate-endpoint-normal.final.spvasm'
    assembly.write_text(transformed)
    recipes=[ [P['tools']['assembler'],'--target-env','spv1.5',str(assembly),'-o',str(module)],
              [P['tools']['validator'],'--target-env','vulkan1.2',str(module)],
              [P['tools']['disassembler'],str(module),'-o',str(readback)] ]
    for argv in recipes:subprocess.run(argv,cwd=R,check=True)
    final=readback.read_text();stock=Path(proof['stock_assembly']['resolved_path']).read_text()
    assert h['interface_signature'](stock)==h['interface_signature'](final)
    modes=lambda text:sorted(re.findall(r'OpExecutionMode %\S+ (.*)$',text,re.M))
    assert modes(stock)==modes(final)
    def barriers(text):
        defs=dict(re.findall(r'^\s*(%\S+) = (.+)$',text,re.M))
        return [[op,*[defs[word] for word in operands.split()]] for op,operands in re.findall(r'^\s*(OpControlBarrier|OpMemoryBarrier) (.+)$',text,re.M)]
    assert barriers(stock)==barriers(final)
    assert sorted(re.findall(r'OpCapability (\S+)',final))==['RoundingModeRTE','Shader']
    assert set(re.findall(r'OpTypeInt (\d+) [01]',final))=={'32'} and re.findall(r'OpTypeFloat (\d+)',final)==['32']
    assert re.search(r'OpDecorate %proposal NoContraction(?:\s|$)',final)
    assert len(re.findall(r'= OpFDiv %float ',final))==1
    assert not re.search(r'\bOpConvertFTo[US]\b',final)
    static()
    result={'schema':'sirius-fp32-estimate-policy-recovery-result-v1','passed':True,'completed_stop':proof,
            'exact_annotation_added':'OpDecorate %proposal NoContraction','one_annotation_and_original_two_RTE_lines_reversal_equal_raw':True,
            'unchanged_frozen_private_and_original_normal_policies_passed':True,'shared_bound_bytes':12596,
            'complete_type_interface_execution_modes_and_ordered_barriers_equal_stock':True,'recipes':recipes,
            'outputs':{str(p):record(p) for p in (assembly,module,readback)},
            'old_failed_passed':False,'old_failed_child_exit':1,'compiler_or_SSA_or_old_disassembly_repeated':False,
            'structural_graph_accepted':False,'runtime_or_numerical_admission':False,'benefit_established':False}
    write(W/'policy-result.json',result)
    write(W/'finish-policy/output-seal.json',{'label':'finish-policy','source_freeze':record(W/'source-freeze.json'),'operator_context':record(W/'operator-context.json'),
          'inputs_before':record(W/'finish-policy/inputs-before.json'),'result':record(W/'policy-result.json'),
          'outputs':{str(p):record(p) for p in (assembly,module,readback)},'structural_graph_accepted':False})
    print('Narrow policy annotation/assembly/validation complete; original failed owner unchanged.',flush=True)

def review():
    static();proof=completed();folder=W/'finish-policy';owner=json.loads((folder/'owner.json').read_text())
    assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['observed_births'] and all(b['absent'] for b in owner['observed_births'])
    assert owner['source']=={'revision':P['operator_revision'],'status':''} and owner['limits']==P['limits']
    action=json.loads((W/'actions.json').read_text())['finish-policy'];assert owner['argv']==action['argv']
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    for name,pin in json.loads((folder/'inputs-before.json').read_text()).items():
        p=Path(name);p=R/p if not p.is_absolute() else p
        assert content(p)=={k:pin[k] for k in ('bytes','sha256')} and str(p.resolve(strict=True))==pin['resolved_path']
    seal=json.loads((folder/'output-seal.json').read_text())
    for key,p in [('source_freeze',W/'source-freeze.json'),('operator_context',W/'operator-context.json'),('inputs_before',folder/'inputs-before.json'),('result',W/'policy-result.json')]:assert seal[key]==record(p)
    expected={str(B/name):record(B/name) for name in ('candidate-endpoint-normal.spvasm','candidate-endpoint-normal.spv','candidate-endpoint-normal.final.spvasm')}
    assert seal['outputs']==expected and seal['structural_graph_accepted'] is False
    result=json.loads((W/'policy-result.json').read_text());assert result['passed'] and result['completed_stop']==proof and result['outputs']==expected
    assert result['old_failed_passed'] is False and result['old_failed_child_exit']==1 and not result['structural_graph_accepted'] and not result['runtime_or_numerical_admission']
    write(W/'policy-acceptance.json',{'schema':'sirius-fp32-estimate-policy-recovery-acceptance-v1','passed':True,'reader_owner':record(folder/'owner.json'),'output_seal':record(folder/'output-seal.json'),'policy_result':record(W/'policy-result.json'),'original_failed_owner':proof['old_failed_owner'],'original_failed_passed':False,'original_failed_child_exit':1,'remaining_assembly_validation_complete':True,'structural_graph_accepted':False,'runtime_or_numerical_admission':False})
    print('Terminal narrow policy recovery accepted; actual graph and lowering admission remain pending.',flush=True)

if __name__=='__main__':
    {'pin-context':pin_context,'finish':finish,'review':review}[sys.argv[1]]()
