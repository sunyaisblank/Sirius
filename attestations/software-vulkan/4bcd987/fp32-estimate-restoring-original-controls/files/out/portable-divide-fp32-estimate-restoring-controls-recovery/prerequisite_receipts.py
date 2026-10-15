"""Selected immutable restoring arithmetic/profile and canonical producer joins.

Historical source/input maps are archive bytes, never replayed as retired paths.
No compiler, parser, policy, tool, provider or test is launched here.
"""
from pathlib import Path
import hashlib,json,os,re,stat

def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'resolved_path':str(path.resolve(strict=True))}

def content(path):return {k:v for k,v in record(path).items() if k!='resolved_path'}

def accepted(folder,revision):
    folder=Path(folder);owner=json.loads((folder/'owner.json').read_text())
    assert owner['status']=='terminal' and owner['passed'] is True and owner['child_exit']==0
    assert owner['source']=={'revision':revision,'status':''} and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    return owner

def archived(pin,path,original):
    path=Path(path);assert path.is_file() and not path.is_symlink()
    assert content(path)=={k:pin[k] for k in ('bytes','sha256')},str(path)
    if 'resolved_path' in pin:assert pin['resolved_path']==str(original)
    if 'mode' in pin:assert stat.S_IMODE(path.stat().st_mode)==pin['mode']
    if 'symlink_target' in pin:assert pin['symlink_target'] is None

def archive(row):
    assert re.fullmatch('[0-9a-f]{40}',row['commit']), 'immutable archive pending'
    root=Path(row['path']);assert content(root/'manifest.json')==row['manifest']
    files=json.loads((root/'manifest.json').read_text())['files']
    for name,pin in row['selected_members'].items():
        path=root/name;assert not path.is_symlink() and content(path)==files[name]==pin,name
    return record(root/'manifest.json')

def verify_canonical_producer(P,R):
    manifests={key:archive(P['archives'][key]) for key in ('apparatus','canonical','arithmetic')}
    arith=Path(P['arithmetic_work']);build=json.loads((arith/'build-review.json').read_text())
    proof=build['canonical_numeric_object_relink']
    assert build['passed'] and build['private_Int64_or_Float64_device_object_used'] is False
    assert record(P['device_object'])==proof['canonical_device_object']
    assert record(P['base_archive'])==proof['canonical_base']
    prior=Path(P['canonical_link_work']);accepted(prior,P['canonical_link_revision'])
    assert record(prior/'owner.json')==proof['canonical_prior_link_owner']
    before=json.loads((prior/'inputs-before.json').read_text())
    for row in P['canonical_producer_rows']:
        pin=before[row['original_input_key']]
        assert content(row['actual'])==content(row['canonical_current'])=={k:pin[k] for k in ('bytes','sha256')}
    current=json.loads((Path(P['apparatus_work'])/'link-candidate/inputs-before.json').read_text())
    accepted(Path(P['apparatus_work'])/'link-candidate',P['runner_preparation_revision'])
    for row in P['matched_link_inputs']:
        assert current[row['original_input_key']]==row['original_pin']
        assert content(row['path'])=={k:row['original_pin'][k] for k in ('bytes','sha256')}
    # Fresh compilation consumes the same original TUs and non-generated .d
    # dependencies, with only source/owned-header namespace relocation.
    for label in ('compile-compute-sealed','compile-observer'):
        accepted(Path(P['apparatus_work'])/label,P['runner_preparation_revision'])
    for row in P['matched_compile_dependencies']:
        before=json.loads(Path(row['receipt']).read_text());assert before[row['original_input_key']]==row['original_pin']
        assert content(row['actual'])=={k:row['original_pin'][k] for k in ('bytes','sha256')}
    return {'device_object':record(P['device_object']),'base_archive':record(P['base_archive']),
            'canonical_prior_link_owner':record(prior/'owner.json'),'accepted_arithmetic_build':record(arith/'build-review.json'),
            'matched_link_input_count':len(P['matched_link_inputs']),'matched_compile_dependency_count':len(P['matched_compile_dependencies']),
            'private_Int64_or_Float64_object_required':False,'fresh_enabled55_argument_capture_claimed':False,
            'manifests':manifests,
            'scope':'Exact canonical Vulkan object/base and unchanged integrated test objects/libraries/source dependencies. Canonical feature enablement remains source/accepted producer evidence; new numerical owners must observe actual three provider maps. No fresh55-field argument equivalence claim.'}

def prior_precontext_stop(P,R):
    row=P['archives']['precontext_stop'];manifest=archive(row)
    work=Path(P['precontext_stop_work']);old=R/'out/portable-divide-fp32-estimate-restoring-controls'
    first=json.loads((work/'root-first-stop.json').read_text())
    independent=json.loads((work/'independent-stop-review.json').read_text())
    root=json.loads((work/'root-context-stop-admission.json').read_text())
    assert first['schema']=='sirius-fp32-restoring-controls-root-precontext-stop-v1'
    assert root['schema']=='sirius-fp32-restoring-controls-root-context-stop-admission-v1'
    assert first['passed'] is True and independent['passed'] is True and root['passed'] is True
    assert first['source_revision']==independent['source_revision']==root['source_revision']==P['runner_preparation_revision']
    assert first['classification']==root['classification']==independent['actual_stop']['classification']=='PREREQUISITE_SCHEMA_MISMATCH_PRECONTEXT_STOP'
    for value in (first,root,independent['actual_stop']):
        assert value['context_attempt_passed'] is False and value['operator_context_accepted'] is False
    assert first['tool_exit_code']==1 and first['tool_session_terminal'] is True
    assert first['producer_compiler_provider_control_or_reader_workloads_launched']==0
    assert first['prospective_actions_unconsumed']==15 and first['prospective_commands_unexecuted_after_first']==17
    assert first['frozen_source_current'] is True and first['source_freeze_count_current']==1095
    assert first['all_absent'] and all(value is True for value in first['all_absent'].values())
    assert first['actual_arithmetic_schema']=='sirius-fp32-estimate-restoring-original-arithmetic-result-v1'
    assert first['expected_schema_literal']=='sirius-fp32-restoring-original-arithmetic-result-v1'
    assert first['accepted_arithmetic_passed'] is True
    assert record(Path(P['arithmetic_work'])/'arithmetic-result.json')==first['accepted_arithmetic_result']
    for key,name in [('root_first_stop','root-first-stop.json'),('independent_stop_review','independent-stop-review.json')]:
        archived(root[key],work/name,old/name)
    for key,name in [('source_freeze','source-freeze.json'),('actions','actions.json'),('packet','packet.json'),('preparation','preparation.json')]:
        archived(first[key],work/name,old/name)
    for pin in first['historical_source_admissions']:
        original=Path(pin['resolved_path']);assert original.parent==old
        archived(pin,work/original.name,original)
        assert json.loads((work/original.name).read_text())['passed'] is True
    # Bind the exact immediately consumed source/preimages, not old1095 current
    # paths. The original SOURCE reports stay unchanged despite the actual stop.
    for name in ('gate.py','prerequisite_receipts.py','verify_producer.py','packet.json','preparation.json','source-freeze.json','actions.unpinned.json','provenance.json'):
        assert (Path(P['recovery_work'])/'originals'/('before-precontext-stop-'+name)).read_bytes()==(work/name).read_bytes()
    return {'manifest':manifest,'root_first_stop':record(work/'root-first-stop.json'),
            'independent_stop_review':record(work/'independent-stop-review.json'),
            'root_stop_admission':record(work/'root-context-stop-admission.json'),
            'context_attempt_passed':False,'operator_context_accepted':False,'tool_exit_code':1,
            'source_admissions_preserved':True,'old1095_paths_replayed_as_current':False,
            'scope':'Accepted stop evidence only. Actual pin-context refused before context/owner/build/provider/test. Separate corrected schema/new namespace changes no original numerical or governor gate.'}

def verify(P,R):
    prior_stop=prior_precontext_stop(P,R)
    manifests={key:archive(row) for key,row in P['archives'].items()}
    for row in P['reader_reference_aliases']:
        assert content(row['original'])==content(row['backed'])
    arith=Path(P['arithmetic_work']);load=lambda n:json.loads((arith/n).read_text())
    ar=load('arithmetic-result.json')
    assert ar['passed'] and ar['schema']=='sirius-fp32-estimate-restoring-original-arithmetic-result-v1'
    assert ar['candidate_macro']==P['candidate_macro']=='SIRIUS_PB32_DIVIDE_FP32_ESTIMATE=1'
    assert ar['original_records']==174268 and ar['original_operation_count']==14 and ar['complete_output_words_per_provider']==697072
    assert ar['fresh_candidate_host_and_full_gpu_plus_reused_stock_equal_original_expected'] and ar['stock_result_reused']
    for n in ['independent-terminal-review.json','root-terminal-admission.json']:assert load(n)['passed'] and load(n)['source_revision']==P['runner_preparation_revision']
    for label in ('execute-host','execute-gpu-candidate'):accepted(arith/label,P['runner_preparation_revision'])
    old_arith=R/'out/portable-divide-fp32-estimate-restoring-arithmetic'
    archived(ar['build_review'],arith/'build-review.json',old_arith/'build-review.json')
    archived(ar['runtime_context'],arith/'runtime-context.json',old_arith/'runtime-context.json')
    arithmetic={'manifest':manifests['arithmetic'],'result':record(arith/'arithmetic-result.json'),
                'terminal_review':record(arith/'independent-terminal-review.json'),'root_admission':record(arith/'root-terminal-admission.json'),
                'scope':'Accepted original174268/14-op/697072-word primitive result and exact stock reuse only; relevant immutable receipts reused, no recursive737 historical current-path replay or retained result inferred.'}
    work=Path(P['profile_work']);original_work=R/'out/portable-divide-fp32-estimate-restoring-profile-reader-recovery'
    result=json.loads((work/'profile-result.json').read_text());approval=json.loads((work/'profile-acceptance.json').read_text())
    assert result['passed'] and approval['passed'] and approval['schema']=='sirius-fp32-restoring-linear-memory-full-profile-acceptance-v1'
    assert approval['original_reader_passed'] is False and approval['original_reader_child_exit']==-15 and approval['original_reader_promoted'] is False
    assert approval['all35_modules_40main_2owner_arrays_7tables_accepted'] and approval['named_graphs']==13 and approval['inlined_camera_footprints']==4
    assert approval['new_disassembly_count']==16 and approval['reused_disassembly_count']==54
    assert approval['arithmetic_or_integrated_or_benefit_or_adoption_admitted'] is False
    assert result['original_failed_reader_promoted'] is False
    for n in ['independent-terminal-review.json','root-terminal-admission.json']:
        v=json.loads((work/n).read_text());assert v['passed'] and v['source_revision']==P['runner_preparation_revision']
    for key,relative in [('profile_result','profile-result.json'),('reader_owner','read-full-profile/owner.json'),('output_seal','read-full-profile/output-seal.json')]:
        archived(approval[key],work/relative,original_work/relative)
    owner=accepted(work/'read-full-profile',P['runner_preparation_revision'])
    assert owner['limits']==P['limits']
    seal=json.loads((work/'read-full-profile/output-seal.json').read_text())
    for key,relative in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before','read-full-profile/inputs-before.json')]:
        archived(seal[key],work/relative,original_work/relative)
    for name,pin in seal['outputs'].items():
        original=Path(name);assert original.is_relative_to(original_work)
        archived(pin,work/original.relative_to(original_work),original)
    profile=result['profile'];assert profile['passed'] and profile['changed_modules']==17 and profile['unchanged_modules']==18
    assert profile['main_arrays_changed']==17 and profile['main_arrays_unchanged']==23 and profile['owner_arrays_unchanged']==2
    assert profile['all_seven_tables_and_nonarray_ABI_constants_equal'] and len(profile['modules'])==35
    assert sorted(row['module'] for row in profile['modules'] if row['changed'])==P['changed_modules']
    assert set(profile['headers']['retained_kernels.h']['changed_arrays'])==set(P['changed_arrays'])
    products=Path(P['candidate_products']).parent
    for row in profile['modules']:
        for mode in ('stock','candidate'):assert record(products/mode/row['module'])==row[mode]
    for name,row in profile['headers'].items():
        for mode in ('stock','candidate'):assert record(products/mode/name)==row[mode]
    chain=result['accepted_compiler_chain'];assert chain==approval['immutable_compiler_chain']
    assert chain['accepted_candidate_products']==90 and chain['accepted_stock_products']==39
    assert chain['original_reader_passed'] is False and chain['original_reader_child_exit']==-15
    assert chain['original_reader_stop_reason']=='diagnostic_sampled_rss_guard'
    assert chain['historical467_source_and1113_maps_replayed_as_current'] is False
    compiled=Path(P['compiled_work']);old_compiled=R/'out/portable-divide-fp32-estimate-restoring-profile-recovery'
    for label,rows in chain['completed_owners'].items():
        accepted(compiled/label,P['runner_preparation_revision'])
        for key,relative in [('owner','owner.json'),('inputs_before','inputs-before.json'),('inputs_after','inputs-after.json'),('output_seal','output-seal.json')]:
            assert record(compiled/label/relative)=={k:rows[key][k] for k in ('bytes','sha256','resolved_path')}
    failed=json.loads((compiled/'read-profile/owner.json').read_text())
    assert failed['passed'] is False and failed['child_exit']==-15 and failed['stop_reason']=='diagnostic_sampled_rss_guard'
    assert content(compiled/'read-profile/owner.json')=={k:chain['original_failed_reader']['owner'][k] for k in ('bytes','sha256')}
    for label in ('reuse-stock','compile-candidate'):
        oldseal=json.loads((compiled/label/'output-seal.json').read_text())
        for name,pin in oldseal['outputs'].items():
            original=Path(name)
            if original.is_relative_to(P['old_products']):actual=products/original.relative_to(P['old_products'])
            else:
                assert original.is_relative_to(old_compiled);actual=compiled/original.relative_to(old_compiled)
            archived(pin,actual,original)
    oldheader=compiled/'candidate/src/sirius/kernels/portable_binary32.h'
    assert content(oldheader)=={k:ar['candidate_header'][k] for k in ('bytes','sha256')}
    source_deltas=[]
    for name,row in P['source_base_rows'].items():
        original=Path(row['stock']);private=Path(row['private']);current=R/'src/sirius/kernels'/name
        assert content(current)==row['current_pin']
        if original.read_bytes()!=current.read_bytes():
            source_deltas.append(name)
            known=P['known_outside_retained_source_deltas'][name]
            assert content(current)==known['canonical4498'] and content(private)==known['original_private']
        if name!='portable_binary32.h':assert original.read_bytes()==private.read_bytes()
        else:assert content(private)==content(oldheader)
        if name in P['retained_source_closure']:assert current.read_bytes()==original.read_bytes()
    assert sorted(source_deltas)==sorted(P['known_outside_retained_source_deltas']) and len(source_deltas)==4 and len(P['source_base_rows'])==33 and len(P['retained_source_closure'])==14
    return {'prior_precontext_stop':prior_stop,'arithmetic':arithmetic,'full_profile_manifest':manifests['profile'],'full_profile_acceptance':record(work/'profile-acceptance.json'),
            'profile_result':record(work/'profile-result.json'),'profile_terminal_review':record(work/'independent-terminal-review.json'),
            'profile_root_admission':record(work/'root-terminal-admission.json'),'compiled_manifest':manifests['compiled'],
            'original_failed_reader':record(compiled/'read-profile/owner.json'),'original_reader_passed':False,'original_reader_child_exit':-15,
            'direct_candidate_products':str(P['candidate_products']),'source_base':'Original33 snapshots/14 actual retained closure; exact4 current source deltas remain outside retained closure.',
            'topology':'2historicalpassedproducer receipts+1separateacceptedprofile; oldRSSreaderfalse/-15, then5fresh producers+8originalcontrols+2readers.',
            'scope':'Exact immutable restoring primitive/profile prerequisites. No original8, benefit, physical-domain, fullscience, rendering or release qualification is inferred.'}
