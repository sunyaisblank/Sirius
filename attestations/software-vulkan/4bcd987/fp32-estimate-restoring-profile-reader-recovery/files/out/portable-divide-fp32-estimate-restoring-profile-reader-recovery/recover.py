"""One root-owned full-profile finisher over immutable accepted compiler bytes.

No stock staging, compiler, SSA, policy, assembler, validator, provider or test is
run. Only missing disassembly and the unchanged profile predicates execute in
the prospective owned reader. Original RSS failure remains false.
"""
from pathlib import Path
import hashlib,json,os,resource,stat,subprocess,sys
sys.dont_write_bytecode=True
W=Path(__file__).resolve().parent
R=W.parents[1]
P=json.loads((W/'packet.json').read_text())


def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'mode':stat.S_IMODE(path.stat().st_mode),
            'resolved_path':str(path.resolve(strict=True)),
            'symlink_target':os.readlink(path) if path.is_symlink() else None}


def content(path):
    actual=record(path)
    return {key:actual[key] for key in ('bytes','sha256')}


def write(path,value):
    path=Path(path);assert not path.exists(),path
    path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')


def static():
    for name,pin in json.loads((W/'source-freeze.json').read_text())['pins'].items():assert record(name)==pin,name


def validate_actions():
    actions=json.loads((W/'actions.json').read_text())
    assert set(actions)=={'read-full-profile'}
    action=actions['read-full-profile']
    assert set(action)=={'kind','revision','argv','frozen_inputs','limits','scope','expected_exit','software_provider'}
    assert action['revision']==P['revision'] and action['kind']=='build'
    assert action['argv']==['/usr/bin/python3','-B',str(W/'recover.py'),'read-profile']
    assert type(action['expected_exit']) is int and action['expected_exit']==0
    assert action['software_provider'] is False
    assert isinstance(action['scope'],str) and action['scope']
    assert action['limits']==P['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
    assert all(type(x) is int for x in action['limits'].values())
    files=action['frozen_inputs']
    assert isinstance(files,list) and len(files)==len(set(files))
    assert all(isinstance(name,str) and Path(name).is_absolute() for name in files)
    assert {str(W/name) for name in ('actions.json','source-freeze.json','operator-context.json')}<=set(files)
    assert not set(files)&set(P['outputs'])
    assert (W/'actions.json').read_bytes()==(W/'actions.unpinned.json').read_bytes()
    return action


def context():
    validate_actions();saved=json.loads((W/'operator-context.json').read_text())
    assert saved['revision']==P['revision']
    for key,name in [('source_freeze','source-freeze.json'),('actions','actions.json'),('preparation','preparation.json')]:assert saved[key]==record(W/name)
    return saved


def completed():
    archive=Path(P['compiled_archive']);old=archive/P['old_work_member'];products=Path(P['original_products'])
    assert content(archive/'manifest.json')==P['compiled_manifest']
    manifest=json.loads((archive/'manifest.json').read_text())['files']
    for member,pin in P['selected_members'].items():
        path=archive/member
        assert not path.is_symlink() and content(path)==manifest[member]==pin,member
    load=lambda name:json.loads((old/name).read_text())
    packet=load('packet.json');first=load('root-first-stop.json')
    independent=load('independent-stop-review.json');root=load('root-reader-rss-stop-admission.json')
    assert first['passed'] and independent['passed'] and root['passed']
    assert first['classification']==root['classification']=='PROFILE_READER_RSS_GUARD_STOP'
    assert first['source_revision']==root['source_revision']==P['revision']
    assert root['original_reader_passed'] is False and root['original_reader_child_exit']==-15
    assert root['original_reader_wrapper_exit']==1 and root['completed_compiler_accepted']
    assert load('root-source-admission.json')['passed'] and load('independent-source-review.json')['passed']
    oldwork=Path(P['old_work']);oldproducts=Path(P['old_products'])
    def historical(expected,path,original):
        assert not path.is_symlink() and content(path)=={key:expected[key] for key in ('bytes','sha256')}
        assert expected['resolved_path']==str(original)
        if 'symlink_target' in expected:assert expected['symlink_target'] is None
        if 'mode' in expected:assert record(path)['mode']==expected['mode']
    historical(root['root_first_stop'],old/'root-first-stop.json',oldwork/'root-first-stop.json')
    historical(root['independent_stop_review'],old/'independent-stop-review.json',oldwork/'independent-stop-review.json')
    historical(root['reader_owner'],old/'read-profile/owner.json',oldwork/'read-profile/owner.json')
    historical(root['compiler_acceptance'],old/'compiler-producer-acceptance.json',oldwork/'compiler-producer-acceptance.json')
    frozen=load('source-freeze.json')['pins'];assert len(frozen)==467
    assert content(old/'source-freeze.json')==P['old_source_freeze']
    receipts={}
    for label in ('reuse-stock','compile-candidate','read-profile'):
        owner=load(label+'/owner.json');failed=label=='read-profile'
        assert owner['status']=='terminal' and owner['passed'] is (not failed) and owner['child_exit']==(-15 if failed else 0)
        assert owner['source']=={'revision':P['revision'],'status':''} and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
        assert owner['stop_reason']==('diagnostic_sampled_rss_guard' if failed else None)
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
        assert owner['limits']==P['limits']
        assert (old/label/'inputs-before.json').read_bytes()==(old/label/'inputs-after.json').read_bytes()
        before=load(label+'/inputs-before.json');assert len(before)==P['old_input_counts'][label]
        receipts[label]={'owner':record(old/label/'owner.json'),'inputs_before':record(old/label/'inputs-before.json'),
                         'inputs_after':record(old/label/'inputs-after.json')}
        if failed:
            assert owner['peak_sampled_rss_kib']>P['limits']['rss_kib']
            assert not (old/label/'output-seal.json').exists()
            continue
        seal=load(label+'/output-seal.json')
        assert seal['label']==label and seal['recipe']==packet['recipes'][label]
        for key,relative in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before',label+'/inputs-before.json')]:
            historical(seal[key],old/relative,oldwork/relative)
        assert set(seal['outputs'])==set(packet['outputs'][label])
        for name,pin in seal['outputs'].items():
            original=Path(name)
            if original.is_relative_to(oldproducts):actual=products/original.relative_to(oldproducts)
            else:
                assert original.is_relative_to(oldwork);actual=old/original.relative_to(oldwork)
            historical(pin,actual,original)
        receipts[label]['output_seal']=record(old/label/'output-seal.json')
    compiler=load('compiler-products.json');acceptance=load('compiler-producer-acceptance.json')
    assert compiler['passed'] and acceptance['passed'] and compiler['compiler_policy_and_validation_completed']
    assert compiler['structural_graph_accepted'] is False and acceptance['structural_graph_accepted'] is False
    assert acceptance['provider_or_numerical_or_benefit_executed'] is False
    historical(acceptance['compiler_result'],old/'compiler-products.json',oldwork/'compiler-products.json')
    for key,relative in [('owner','compile-candidate/owner.json'),('inputs_before','compile-candidate/inputs-before.json'),
                         ('inputs_after','compile-candidate/inputs-after.json'),('output_seal','compile-candidate/output-seal.json')]:
        historical(acceptance['compiler_receipt'][key],old/relative,oldwork/relative)
    assert len(compiler['products'])==90 and len(list((products/'candidate').iterdir()))==90
    assert len(list((products/'stock').iterdir()))==39
    compile_before=load('compile-candidate/inputs-before.json')
    for relative in P['compiled_source_members']:
        original=oldwork/relative;key=str(original.relative_to(R));path=old/relative
        historical(compile_before[key],path,original)
        assert compile_before[key]=={field:frozen[str(original)][field] for field in compile_before[key]}
    assert content(old/'candidate/src/sirius/kernels/portable_binary32.h')==P['candidate_header']
    # The old reader actually consumed these exact stock/candidate bytes. Join
    # only relevant sealed products/receipts, not all1113 historical current paths.
    before=load('read-profile/inputs-before.json')
    for mode in ('stock','candidate'):
        for path in (products/mode).iterdir():
            original=oldproducts/mode/path.name;key=str(original.relative_to(R))
            historical(before[key],path,original)
    partial=archive/P['partial_assembly_member']
    assert content(partial)==P['partial_assembly_pin']
    historical(first['partial_assembly'],partial,oldwork/'analysis/candidate/retained_camera_portable.spvasm')
    assert first['partial_assembly_accepted_as_profile'] is False
    assert not (old/'profile-result.json').exists() and not (old/'profile-acceptance.json').exists()
    assert P['changed_modules']==packet['changed_modules'] and P['changed_arrays']==packet['changed_arrays']
    assert P['module_arrays']==packet['module_arrays']
    assert set(P['saved_assembly'])==set(packet['saved_assembly'])|{'candidate/retained_camera_portable.spv'}
    for key,row in P['saved_assembly'].items():
        mode,name=key.split('/',1)
        assert Path(row['source_module'])==products/mode/name
        assert content(row['source_module'])==row['source_module_pin']
        assert content(row['path'])==row['pin']
        if key in packet['saved_assembly']:
            original=packet['saved_assembly'][key]
            assert row['path']==original['path'] and row['pin']==original['pin']
            assert row['source_module_pin']==original['source_module_pin']
        else:assert Path(row['path'])==partial
    # Preserve the original35 accepted stock ASM→binary producer joins, plus
    # the18 exact unchanged candidate aliases. No old policy/reader is rerun.
    saved_proof=P['saved_original_assembly_proof']
    assert saved_proof==packet['assembly_archive']
    proof_archive=Path(saved_proof['path'])
    assert content(proof_archive/'manifest.json')==saved_proof['manifest']
    proof_manifest=json.loads((proof_archive/'manifest.json').read_text())['files']
    for member,pin in saved_proof['selected_members'].items():
        assert content(proof_archive/member)==proof_manifest[member]==pin
    proof_work=proof_archive/saved_proof['work_member']
    proof_result=json.loads((proof_work/'profile-result.json').read_text())
    proof_acceptance=json.loads((proof_work/'profile-acceptance.json').read_text())
    assert proof_result['passed'] and proof_acceptance['passed']
    assert json.loads((proof_work/'independent-terminal-review.json').read_text())['passed']
    assert {key:proof_acceptance['profile_result'][key] for key in ('bytes','sha256')}==content(proof_work/'profile-result.json')
    proof_rows={row['module']:row for row in proof_result['profile']['modules']};assert len(proof_rows)==35
    for name in P['module_names']:
        original=proof_rows[name]['assembly_binary_joins']['stock']
        row=P['saved_assembly']['stock/'+name]
        assert row['pin']=={key:original['assembly'][key] for key in ('bytes','sha256')}
        assert row['source_module_pin']=={key:original['source_module'][key] for key in ('bytes','sha256')}
        historical_asm=Path(original['assembly']['resolved_path']);actual_asm=Path(row['path'])
        prefix=R/'out/portable-divide-fp64-estimate-profile-recovery'
        assert actual_asm==historical_asm or (historical_asm.is_relative_to(prefix) and actual_asm==proof_work/historical_asm.relative_to(prefix))
    for name in P['unchanged_helpers']:
        assert (W/name).read_bytes()==(old/name).read_bytes(),name
    delta=json.loads((W/'source-delta.json').read_text())
    for name,replacements in delta['exact_reversal'].items():
        restored=(W/name).read_text()
        for row in replacements:
            assert restored.count(row['new'])==1
            restored=restored.replace(row['new'],row['old'],1)
        assert restored==(old/name).read_text(),name
    return {'immutable_manifest':record(archive/'manifest.json'),'accepted_compiler':record(old/'compiler-producer-acceptance.json'),
            'completed_owners':{key:row for key,row in receipts.items() if key!='read-profile'},
            'original_failed_reader':receipts['read-profile'],'original_reader_passed':False,'original_reader_child_exit':-15,
            'original_reader_stop_reason':'diagnostic_sampled_rss_guard','root_stop_admission':record(old/'root-reader-rss-stop-admission.json'),
            'independent_stop_review':record(old/'independent-stop-review.json'),'prior_source_freeze':record(old/'source-freeze.json'),
            'accepted_candidate_products':90,'accepted_stock_products':39,
            'partial_dump':record(partial),'partial_dump_prior_graph_acceptance':False,
            'historical467_source_and1113_maps_replayed_as_current':False,
            'scope':'Exact selected immutable accepted compiler/seals and original failed RSS reader. Saved policy/compiler/source admissions reused without tool or current-path replay; partial dump must pass fresh whole-module alpha verification.'}


def pin_context():
    validate_actions();static();completed()
    assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=R)
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()==P['revision']
    assert not (W/'operator-context.json').exists() and not (W/'read-full-profile').exists() and not (W/'analysis').exists()
    assert all(not Path(name).exists() for name in P['outputs']) and not (W/'profile-acceptance.json').exists()
    write(W/'operator-context.json',{'revision':P['revision'],'source_freeze':record(W/'source-freeze.json'),
          'actions':record(W/'actions.json'),'preparation':record(W/'preparation.json')})


def read_profile():
    validate_actions();static();context();proof=completed()
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));_,prior=resource.getrlimit(resource.RLIMIT_FSIZE)
    cap=64*1024**2;limit=cap if prior==resource.RLIM_INFINITY else min(prior,cap)
    resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
    assert all(not Path(name).exists() for name in P['outputs'])
    from products_profile import products_review
    result=products_review(P,R,W,Path(P['original_products']),Path(P['original_work']))
    assert result['changed_modules']==17 and result['unchanged_modules']==18
    generated={str(path):record(path) for path in (W/'analysis/candidate').glob('*.spvasm')}
    assert set(generated)==set(P['new_assembly_outputs']) and len(generated)==16
    partial=result['modules'][P['module_names'].index('retained_camera_portable.spv')]
    assert partial['assembly_binary_joins']['candidate']['reused_saved_assembly']
    assert partial['divide_footprint']['persisted_compiler_policy_join']['complete_assembler_readback_instructions_equal']>0
    write(W/'profile-result.json',{'schema':'sirius-fp32-restoring-linear-memory-full-profile-result-v1','passed':True,
          'profile':result,'accepted_compiler_chain':proof,'new16_assembly':generated,'reused54_assembly':P['saved_assembly'],
          'dominance_relation':'Entry-reachable target becomes unreachable when the queried dominating block is removed; boolean-only cached answers.',
          'original_failed_reader_promoted':False,'arithmetic_or_integrated_or_benefit_or_adoption_admitted':False})
    static()
    write(W/'read-full-profile/output-seal.json',{'label':'read-full-profile','source_freeze':record(W/'source-freeze.json'),
          'operator_context':record(W/'operator-context.json'),'inputs_before':record(W/'read-full-profile/inputs-before.json'),
          'outputs':{name:record(name) for name in P['outputs']}})


def review():
    validate_actions();static();context();proof=completed();folder=W/'read-full-profile'
    action=json.loads((W/'actions.json').read_text())['read-full-profile'];owner=json.loads((folder/'owner.json').read_text())
    assert owner['status']=='terminal' and owner['passed'] and owner['child_exit']==0
    assert owner['source']=={'revision':P['revision'],'status':''} and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['argv']==action['argv'] and owner['limits']==P['limits']
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    before=json.loads((folder/'inputs-before.json').read_text())
    for name in action['frozen_inputs']:
        path=Path(name);key=str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
        assert before[key]=={k:record(path)[k] for k in before[key]},name
    seal=json.loads((folder/'output-seal.json').read_text())
    for key,path in [('source_freeze',W/'source-freeze.json'),('operator_context',W/'operator-context.json'),('inputs_before',folder/'inputs-before.json')]:assert seal[key]==record(path)
    assert seal['outputs']=={name:record(name) for name in P['outputs']}
    result=json.loads((W/'profile-result.json').read_text())
    assert result['passed'] and result['accepted_compiler_chain']==proof
    assert result['profile']['changed_modules']==17 and result['profile']['unchanged_modules']==18
    assert result['original_failed_reader_promoted'] is False
    write(W/'profile-acceptance.json',{'schema':'sirius-fp32-restoring-linear-memory-full-profile-acceptance-v1','passed':True,
          'profile_result':record(W/'profile-result.json'),'reader_owner':record(folder/'owner.json'),'output_seal':record(folder/'output-seal.json'),
          'immutable_compiler_chain':proof,'all35_modules_40main_2owner_arrays_7tables_accepted':True,
          'named_graphs':13,'inlined_camera_footprints':4,'new_disassembly_count':16,'reused_disassembly_count':54,
          'original_reader_passed':False,'original_reader_child_exit':-15,'original_reader_promoted':False,
          'arithmetic_or_integrated_or_benefit_or_adoption_admitted':False,
          'camera_scope':'Actual typed five-step guarded dataflow, closed Phi/control and linear-memory exact dominance queries; full inlined normalization/Round/fallback also relies on frozen source/original O1, unchanged from the preceding reader contract.'})


if __name__=='__main__':{'pin-context':pin_context,'read-profile':read_profile,'review':review}[sys.argv[1]]()
