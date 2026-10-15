"""Immutable arithmetic/compiler/profile receipt joins; no tools or provider."""
from pathlib import Path
import hashlib,json,re

def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return {'bytes':path.stat().st_size,'sha256':digest,'resolved_path':str(path.resolve(strict=True))}

def content(path):return {k:v for k,v in record(path).items() if k!='resolved_path'}

def accepted(folder):
    folder=Path(folder);value=json.loads((folder/'owner.json').read_text())
    assert value['status']=='terminal' and value['passed'] and value['child_exit']==0
    assert value['stop_reason'] is None and value['cleanup_errors']==[] and value['remaining_owned_processes']==[]
    assert value['source_unchanged'] and value['whole_inputs_unchanged'] and all(row['absent'] for row in value['observed_births'])
    assert value['observed_births'] and (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    return value

def verify(P,R):
    assert re.fullmatch('[0-9a-f]{40}',P['profile_archive_commit']), 'actual profile archive pending'
    assert re.fullmatch('[0-9a-f]{64}',P['profile_archive_manifest_sha256']), 'actual profile manifest pending'
    def manifest(archive,digest):
        archive=Path(archive);assert content(archive/'manifest.json')['sha256']==digest
        for name,pin in json.loads((archive/'manifest.json').read_text())['files'].items():
            path=archive/name;assert not path.is_symlink() and content(path)==pin,name
        return record(archive/'manifest.json')
    def join(name,pin,archive):
        original=Path(name)
        if not original.is_absolute():original=R/original
        saved=Path(archive)/'files'/original.relative_to(R) if original.is_relative_to(R) else None
        owned=('out/portable-divide-fp64-estimate-arithmetic','bin/linux-gcc/portable-divide-fp64-estimate-arithmetic',
               'out/portable-divide-fp64-estimate-controls','bin/linux-gcc/portable-divide-fp64-estimate-controls',
               'out/portable-divide-fp64-estimate-profile-recovery')
        if any(original.is_relative_to(R/root) for root in owned):assert saved is not None and saved.is_file(),name
        actual=saved if saved is not None and saved.is_file() else original
        assert content(actual)=={key:pin[key] for key in ('bytes','sha256')},name
        if actual==original:assert record(actual)==pin,name
        else:assert not actual.is_symlink() and pin['resolved_path']==str(original),name
        return actual
    arithmetic=Path(P['arithmetic_archive']);arith=arithmetic/'files/out/portable-divide-fp64-estimate-arithmetic'
    arith_manifest=manifest(arithmetic,P['arithmetic_manifest_sha256'])
    result_path=arith/'arithmetic-result.json';assert content(result_path)['sha256']==P['arithmetic_inner_result_sha256']
    ar=json.loads(result_path.read_text());assert ar['passed'] and ar['schema']=='sirius-fp64-estimate-original-arithmetic-result-v1'
    assert ar['original_records']==174268 and ar['original_operation_count']==14 and ar['complete_output_words_per_provider']==697072
    assert ar['fresh_candidate_host_and_full_gpu_plus_reused_stock_equal_original_expected'] and ar['stock_result_reused']
    assert ar['candidate_macro']==P['candidate_macro']=='SIRIUS_PB32_DIVIDE_FP64_ESTIMATE=1'
    arithmetic_frozen=json.loads((arith/'source-freeze.json').read_text())['pins'];assert len(arithmetic_frozen)==1886
    for name,pin in arithmetic_frozen.items():join(name,pin,arithmetic)
    assert len(ar['owners'])==10
    for label,pin in ar['owners'].items():
        accepted(arith/label);join(str(R/'out/portable-divide-fp64-estimate-arithmetic'/label/'owner.json'),pin,arithmetic)
        for name,pin in json.loads((arith/label/'inputs-before.json').read_text()).items():join(name,pin,arithmetic)
    assert (arith/'execute-gpu-candidate/actual.bin').read_bytes()==(arith/'expected.bin').read_bytes()
    terminal=arith/'independent-terminal-review.json'
    assert content(terminal)['sha256']==P['arithmetic_terminal_review_sha256'] and json.loads(terminal.read_text())['passed']
    arith_receipt={'result':record(result_path),'manifest':arith_manifest,'terminal_review':record(terminal),'all_original_static_and10owner_input_joins':True,'scope':'Original174268/14-op primitive gate with full697072word candidate provider output and exact stock reuse; no retained or hit-rate inference.'}
    compiled=Path(P['compiler_stop_archive']);old=compiled/'files/out/portable-divide-fp64-estimate-controls'
    compiled_manifest=manifest(compiled,P['compiler_stop_manifest_sha256'])
    profile_archive=Path(P['profile_archive']);work=profile_archive/'files/out/portable-divide-fp64-estimate-profile-recovery'
    profile_manifest=manifest(profile_archive,P['profile_archive_manifest_sha256'])
    approval_path=work/'profile-acceptance.json';assert content(approval_path)['sha256']==P['profile_acceptance_sha256']
    approval=json.loads(approval_path.read_text());assert approval['passed'] and approval['schema']=='sirius-fp64-estimate-full-profile-recovery-acceptance-v1'
    assert approval['original_failed_wrapper_passed'] is False and approval['original_failed_wrapper_child_exit']==1
    assert approval['new_reader_owner_passed'] and approval['original35_modules_40main_2owner_arrays_7tables_full_profile_accepted']
    assert approval['new_disassembly_count']==64 and approval['reused_disassembly_count']==6 and not approval['integrated_controls_executed'] and not approval['production_adopted']
    assert approval['original_archive_manifest']==compiled_manifest
    for key,relative in [('profile_result','profile-result.json'),('reader_owner','read-full-profile/owner.json'),('output_seal','read-full-profile/output-seal.json')]:
        assert join(str(R/'out/portable-divide-fp64-estimate-profile-recovery'/relative),approval[key],profile_archive)==work/relative
    accepted(work/'read-full-profile')
    profile_frozen=json.loads((work/'source-freeze.json').read_text())['pins'];assert len(profile_frozen)==2473
    for name,pin in profile_frozen.items():join(name,pin,profile_archive)
    for name,pin in json.loads((work/'read-full-profile/inputs-before.json').read_text()).items():join(name,pin,profile_archive)
    seal=json.loads((work/'read-full-profile/output-seal.json').read_text())
    for key,relative in [('source_freeze','source-freeze.json'),('operator_context','operator-context.json'),('inputs_before','read-full-profile/inputs-before.json'),('profile_result','profile-result.json')]:
        assert join(str(R/'out/portable-divide-fp64-estimate-profile-recovery'/relative),seal[key],profile_archive)==work/relative
    for name,pin in seal['new64_assembly'].items():join(name,pin,profile_archive)
    reviewed=work/'independent-terminal-review.json'
    assert content(reviewed)['sha256']==P['profile_terminal_review_sha256'] and json.loads(reviewed.read_text())['passed']
    result=json.loads((work/'profile-result.json').read_text());assert result['passed']
    profile=result['profile'];assert profile['passed'] and profile['changed_modules']==17 and profile['unchanged_modules']==18
    assert profile['main_arrays_changed']==17 and profile['main_arrays_unchanged']==23 and profile['owner_arrays_unchanged']==2 and profile['all_seven_tables_and_nonarray_ABI_constants_equal']
    assert len(profile['modules'])==35
    assert sorted(row['module'] for row in profile['modules'] if row['changed'])==P['changed_modules']
    assert set(profile['headers']['retained_kernels.h']['changed_arrays'])==set(P['changed_arrays'])
    products=Path(P['candidate_products']).parent
    assert products==compiled/'files/bin/linux-gcc/portable-divide-fp64-estimate-controls'
    for row in profile['modules']:
        for mode in ('stock','candidate'):
            assert record(products/mode/row['module'])==row[mode]
    for name,row in profile['headers'].items():
        for mode in ('stock','candidate'):assert record(products/mode/name)==row[mode]
    proof=result['completed_producers'];assert proof['original_failed_passed'] is False and proof['original_failed_child_exit']==1
    assert len(proof['three_passed_original_producers'])==3 and len(proof['sealed78_products'])==78
    for name,pin in proof['sealed78_products'].items():assert record(name)==pin
    assert approval['original_three_completed_producers']==proof['three_passed_original_producers']
    for label,rows in proof['three_passed_original_producers'].items():
        accepted(old/label)
        for field,relative in [('owner','owner.json'),('inputs_before','inputs-before.json'),('inputs_after','inputs-after.json'),('output_seal','output-seal.json')]:assert record(old/label/relative)==rows[field]
    failed=json.loads((old/'verify-products/owner.json').read_text());assert failed['passed'] is False and failed['child_exit']==1
    assert record(old/'verify-products/owner.json')==proof['original_failed_owner']==approval['original_failed_owner']
    source=json.loads((old/'stock-reuse-acceptance.json').read_text())['retained_source_profile']
    assert source['all14_consumed_source_files_equal_canonical4498'] and len(source['complete_retained_import_include_closure'])==14
    assert source['current_equal_source_count']==29 and source['exact4_current_deltas_outside_retained_closure']==['gr_ellis_topology.slang','gr_integrator.slang','parity_probe.slang','trace.slang']
    for row in source['historical33_original_private_sources'].values():assert record(row['canonical4498']['resolved_path'])==row['canonical4498']
    assert content(old/'candidate/src/sirius/kernels/portable_binary32.h')=={key:ar['candidate_header'][key] for key in ('bytes','sha256')}
    return {'arithmetic':arith_receipt,'compiled_manifest':compiled_manifest,'historical_three_producers':proof['three_passed_original_producers'],
            'full_profile_manifest':profile_manifest,'full_profile_acceptance':record(approval_path),'profile_result':record(work/'profile-result.json'),'profile_terminal_review':record(reviewed),
            'original_failed_owner':proof['original_failed_owner'],'original_failed_wrapper_passed':False,'original_failed_wrapper_child_exit':1,
            'direct_candidate_products':str(Path(P['candidate_products'])),'topology':'3historicalpassedproducers+1separateacceptedprofile, oldfailedreaderfalse/1, then5freshproducer+8originalcontrols+2readers.',
            'scope':'Original full retained profile and arithmetic prerequisite receipts only. Public optionalFloat64 provider compatibility/fallback is not implemented; no integrated/benefit/science qualification inferred.'}
