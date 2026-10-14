"""Root-owned parse-only recovery; original compiler wrapper remains failed."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
sys.dont_write_bytecode=True
W=Path(__file__).resolve().parent
R=W.parents[1]
P=json.loads((W/'packet.json').read_text())

def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'resolved_path':str(path.resolve(strict=True))}

def content(path):return {k:v for k,v in record(path).items() if k!='resolved_path'}

def write(path,value):
    path=Path(path);assert not path.exists(),path
    path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')

def verify_static():
    for name,pin in json.loads((W/'source-freeze.json').read_text())['pins'].items():assert record(name)==pin,name

def terminal(folder,passed=True):
    folder=Path(folder);owner=json.loads((folder/'owner.json').read_text())
    assert owner['status']=='terminal' and owner['passed'] is passed and owner['child_exit']==(0 if passed else 1)
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    return owner

def original_completed():
    archive=Path(P['compiled_archive']);work=Path(P['original_work']);products=Path(P['original_products'])
    assert content(archive/'manifest.json')==P['compiled_manifest']
    manifest=json.loads((archive/'manifest.json').read_text())
    for name,pin in manifest['files'].items():
        path=archive/name;assert not path.is_symlink() and content(path)==pin,name
    def join(name,pin):
        original=Path(name)
        if not original.is_absolute():original=R/original
        saved=archive/'files'/original.relative_to(R) if original.is_relative_to(R) else None
        actual=saved if saved is not None and saved.is_file() else original
        assert content(actual)=={key:pin[key] for key in ('bytes','sha256')},name
        if actual==original:assert record(actual)==pin,name
        else:assert not actual.is_symlink() and pin['resolved_path']==str(original),name
        return actual
    frozen=json.loads((work/'source-freeze.json').read_text())['pins'];assert len(frozen)==2324
    for name,pin in frozen.items():join(name,pin)
    original_packet=json.loads((work/'packet.json').read_text())
    assert P['changed_modules']==original_packet['changed_modules'] and P['changed_arrays']==original_packet['changed_arrays']
    receipts={};sealed_products={}
    for label in ('reuse-products-stock','compile-products-candidate','native-policy-selfcheck'):
        terminal(work/label)
        before=json.loads((work/label/'inputs-before.json').read_text())
        for name,pin in before.items():join(name,pin)
        seal=json.loads((work/label/'output-seal.json').read_text())
        assert seal['label']==label and seal['recipe']==original_packet['recipes'][label]
        for key,relative in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before',label+'/inputs-before.json')]:
            assert join(str(R/'out/portable-divide-fp64-estimate-controls'/relative),seal[key])==work/relative
        assert set(seal['outputs'])==set(original_packet['outputs'][label])
        for name,pin in seal['outputs'].items():
            actual=join(name,pin)
            if actual.is_relative_to(products):sealed_products[str(actual)]=record(actual)
        receipts[label]={'owner':record(work/label/'owner.json'),'inputs_before':record(work/label/'inputs-before.json'),'inputs_after':record(work/label/'inputs-after.json'),'output_seal':record(work/label/'output-seal.json')}
    assert len(sealed_products)==78
    for mode in ('stock','candidate'):
        assert len(list((products/mode).iterdir()))==39
        assert len(list((products/mode).glob('*.spv')))==35
        assert all(str(path) in sealed_products for path in (products/mode).iterdir())
    failed=terminal(work/'verify-products',False)
    assert record(work/'verify-products/owner.json')['sha256']==P['failed_owner_sha256']
    for name,pin in json.loads((work/'verify-products/inputs-before.json').read_text()).items():join(name,pin)
    stderr=(work/'verify-products/stderr.log').read_text()
    assert 'camera_profile.py", line 121, in inspect' in stderr
    assert "assert all(word(k)[1]==u32" in stderr and stderr.rstrip().endswith('AssertionError')
    assert not (work/'verify-products/output-seal.json').exists() and not (work/'products-acceptance.json').exists()
    assert json.loads((work/'independent-terminal-stop-review.json').read_text())['passed']
    assert record(work/'independent-terminal-stop-review.json')['sha256']==P['terminal_stop_review_sha256']
    saved={str(Path(row['path'])) for row in P['saved_assembly'].values()}
    assert saved=={str(path) for path in (products/'analysis').rglob('*.spvasm')} and len(saved)==6
    expected_saved={mode+'/'+name for mode in ('stock','candidate') for name in
                    ('retained_camera.spv','retained_camera_fp64.spv','retained_camera_portable.spv')}
    assert set(P['saved_assembly'])==expected_saved
    for key,row in P['saved_assembly'].items():
        mode,name=key.split('/',1)
        assert Path(row['path'])==products/'analysis'/mode/(name+'asm')
        assert Path(row['source_module'])==products/mode/name
        assert content(row['path'])==row['pin'] and content(row['source_module'])==row['source_module_pin']
    return {'three_passed_original_producers':receipts,'sealed78_products':sealed_products,
            'original_failed_owner':record(work/'verify-products/owner.json'),'original_failed_passed':failed['passed'],'original_failed_child_exit':failed['child_exit'],
            'original_terminal_stop_review':record(work/'independent-terminal-stop-review.json'),
            'scope':'Exact immutable completed stock39/candidate39/source/tool/recipe/input/output/policy reuse. Original failed profile owner remains false/1; no compiler, validator or policy rerun.'}

def signature_controls():
    from camera_profile import integer_division_signature
    stock=Path(P['saved_assembly']['stock/retained_camera_portable.spv']['path']).read_text()
    candidate=Path(P['saved_assembly']['candidate/retained_camera_portable.spv']['path']).read_text()
    expected=P['original_camera_integer_signature']
    assert integer_division_signature(stock)==integer_division_signature(candidate)==expected
    types=dict(re.findall(r'(%\S+) = OpTypeInt 32 ([01])',candidate))
    unsigned=next(k for k,v in types.items() if v=='0');signed=next(k for k,v in types.items() if v=='1')
    def instruction(opcode):
        return next(line for line in candidate.splitlines() if re.fullmatch(r'\s*%\S+ = '+opcode+r' %\S+ %\S+ %\S+',line))
    def once(line,replacement):
        assert candidate.splitlines().count(line)==1
        return candidate.replace(line+'\n',replacement+'\n',1)
    sdiv=instruction('OpSDiv');udiv=instruction('OpUDiv')
    cases={'signed_op_unsigned_result':once(sdiv,sdiv.replace('OpSDiv '+signed+' ','OpSDiv '+unsigned+' ',1)),
           'unsigned_op_signed_result':once(udiv,udiv.replace('OpUDiv '+unsigned+' ','OpUDiv '+signed+' ',1)),
           'wrong_opcode_signed_result':once(sdiv,sdiv.replace('OpSDiv ','OpUDiv ',1)),
           'wide_unsigned_type':candidate.replace(unsigned+' = OpTypeInt 32 0',unsigned+' = OpTypeInt 64 0',1),
           'missing_signed_count':once(sdiv,''),
           'extra_unsigned_count':candidate+'\n%recovery_extra_count = '+udiv.split(' = ',1)[1]+'\n'}
    rows=[]
    for name,text in cases.items():
        assert text!=candidate
        try:
            assert integer_division_signature(text)==expected
        except (AssertionError,KeyError):pass
        else:raise AssertionError('negative signature control accepted: '+name)
        rows.append({'name':name,'bytes':len(text.encode()),'sha256':hashlib.sha256(text.encode()).hexdigest(),'rejected':True})
    return {'passed':True,'positive_complete_modules':2,'negative_controls':rows,'exact_original_signature':expected,
            'scope':'Complete preserved module text, narrow typed-op/count predicate controls. Textual negative mutations are not assembled or claimed valid SPIR-V; no numerical/driver execution.'}

def pin_context():
    verify_static();original_completed()
    assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=R)
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip();assert revision==P['revision']
    assert not (W/'analysis').exists()
    write(W/'operator-context.json',{'revision':revision,'source_freeze':record(W/'source-freeze.json'),'actions':record(W/'actions.json'),'preparation':record(W/'preparation.json')})

def phase():
    verify_static();context=json.loads((W/'operator-context.json').read_text())
    assert context['source_freeze']==record(W/'source-freeze.json') and context['actions']==record(W/'actions.json') and context['preparation']==record(W/'preparation.json')
    proof=original_completed();controls=signature_controls()
    from products_profile import products_review
    profile=products_review(P,R,W,Path(P['original_products']),Path(P['original_work']))
    assert profile['passed']
    generated={str(path):record(path) for path in (W/'analysis').rglob('*.spvasm')}
    assert set(generated)==set(P['new_assembly_outputs']) and len(generated)==64
    write(W/'profile-result.json',{'schema':'sirius-fp64-estimate-full-profile-recovery-result-v1','passed':True,'profile':profile,'completed_producers':proof,'signature_controls':controls,'new64_assembly':generated,'reused6_assembly':P['saved_assembly'],'original_failed_wrapper_passed':False,'original_failed_wrapper_child_exit':1,'integrated_controls_executed':False,'production_adopted':False})
    write(W/'read-full-profile/output-seal.json',{'source_freeze':record(W/'source-freeze.json'),'operator_context':record(W/'operator-context.json'),'inputs_before':record(W/'read-full-profile/inputs-before.json'),'profile_result':record(W/'profile-result.json'),'new64_assembly':generated})

def review():
    verify_static();owner=terminal(W/'read-full-profile');proof=original_completed()
    context=json.loads((W/'operator-context.json').read_text());assert context['source_freeze']==record(W/'source-freeze.json') and context['actions']==record(W/'actions.json') and context['preparation']==record(W/'preparation.json')
    seal=json.loads((W/'read-full-profile/output-seal.json').read_text());assert seal['source_freeze']==record(W/'source-freeze.json') and seal['operator_context']==record(W/'operator-context.json')
    assert seal['inputs_before']==record(W/'read-full-profile/inputs-before.json') and seal['profile_result']==record(W/'profile-result.json')
    action=json.loads((W/'actions.json').read_text())['read-full-profile'];assert owner['argv']==action['argv']
    before=json.loads((W/'read-full-profile/inputs-before.json').read_text())
    for name in action['frozen_inputs']:
        path=Path(name);key=str(path.relative_to(R)) if path.is_relative_to(R) else str(path);assert before[key]==record(path),name
    generated={str(path):record(path) for path in (W/'analysis').rglob('*.spvasm')}
    assert generated==seal['new64_assembly'] and set(generated)==set(P['new_assembly_outputs'])
    result=json.loads((W/'profile-result.json').read_text());assert result['passed'] and result['original_failed_wrapper_passed'] is False and result['original_failed_wrapper_child_exit']==1
    assert result['completed_producers']==proof and result['reused6_assembly']==P['saved_assembly'] and result['new64_assembly']==generated
    assert result['profile']['changed_modules']==17 and result['profile']['unchanged_modules']==18 and result['signature_controls']['passed']
    assert len(result['signature_controls']['negative_controls'])==6 and all(row['rejected'] for row in result['signature_controls']['negative_controls'])
    write(W/'profile-acceptance.json',{'schema':'sirius-fp64-estimate-full-profile-recovery-acceptance-v1','passed':True,'profile_result':record(W/'profile-result.json'),'reader_owner':record(W/'read-full-profile/owner.json'),'output_seal':record(W/'read-full-profile/output-seal.json'),'original_archive_manifest':record(Path(P['compiled_archive'])/'manifest.json'),'original_three_completed_producers':proof['three_passed_original_producers'],'original_failed_owner':proof['original_failed_owner'],'original_failed_wrapper_passed':False,'original_failed_wrapper_child_exit':1,'new_reader_owner_passed':True,'original35_modules_40main_2owner_arrays_7tables_full_profile_accepted':True,'new_disassembly_count':64,'reused_disassembly_count':6,'integrated_controls_executed':False,'production_adopted':False,'scope':'Separate post-terminal parse-only recovery of exact completed compiler outputs. Original failed reader is not promoted. Matched producer/eight original controls require a separate admitted continuation; no numerical/performance/adoption or public-provider compatibility verdict.'})

if __name__=='__main__':{'pin-context':pin_context,'phase':phase,'review':review}[sys.argv[1]]()
