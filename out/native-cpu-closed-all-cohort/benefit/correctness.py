"""Require fresh AX C8 and exact current static composition before any ABBA slot."""
from pathlib import Path
import hashlib
import json


def check():
    root=Path.cwd()
    current=root/'out/native-cpu-closed-all-cohort'
    provider=root/'out/native-cpu-closed-all-provider/provider'
    def doc(path):
        return json.loads(path.read_text())
    def pin(path):
        with path.open('rb') as stream:
            return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())
    manifest_path=root/'attestations/qualification/363aa0c2/native-cpu-closed-all-original-c8/MANIFEST.json'
    assert pin(manifest_path)['sha256']=='7603fbd5ed5f98178550ff78ac04b4b64603aecd70d14cd29a42fb1e2eb397f6'
    manifest=doc(manifest_path)
    assert manifest['status']=='CLOSED_THREE_STAGE_PROVIDER_ORIGINAL_C8_ALL32_PASS_BROADER_GATES_OPEN'
    assert manifest['C8_observed_scalar_stages']==['transport']
    assert manifest['C8_observed_closed_SIMD_stages']==['endpoint','dopri_phase']
    assert manifest['closed_SIMD_stages_unobserved_by_first_packet_C8']==['transport']
    assert len(manifest['owners'])==4
    for name,saved in manifest['owners'].items():
        path=root/name;record=doc(path)
        assert pin(path)==saved['identity'] and record['source']==saved['source']
        assert record['passed'] and record['status']=='terminal' and record['child_exit']==0
        assert record['source_unchanged'] and record['whole_inputs_unchanged']
        assert not record['remaining_owned_processes'] and not record['cleanup_errors']
        assert record['observed_births'] and all(item['absent'] for item in record['observed_births'])
        assert path.with_name('inputs-before.json').read_bytes()==path.with_name('inputs-after.json').read_bytes()
    accepted=doc(current/'control-8-acceptance.json')
    assert pin(current/'control-8-acceptance.json')==manifest['original_C8_acceptance']
    assert accepted['passed'] and accepted['all32_readbacks_match_preserved_originals']
    assert accepted['dispatch_admission']['all_observed_private_rows_have_positive_live_job']
    assert accepted['dispatch_admission']['job_and_descriptor_views_match_actual_packet_arguments']
    optional=accepted['optional_math']
    assert optional['all_selected_optional_math_observed']
    assert optional['genuine_private_FP_scope_observed']
    assert optional['component_numerical_qualification_transferred'] is False
    assert accepted['transport_component_scope']['component_interface']==dict(
        compiled_interface_passed=True,static_symbol_binding_passed=True,
        numerical_qualification_transferred=False,genuine_job_numerical_gate_required=True)
    assert {'endpoint','dopri_phase'}<=set(optional['observed_optional_math_stages'])
    assert set(optional['observed_optional_math_stages'])==set(optional['selected_simd_stages'])
    assert set(optional['selected_simd_stages'])=={'endpoint','dopri_phase'}
    assert accepted['first_raw_SIMD_entries_completed']==2
    assert accepted['observed_Transport12_implementation']=='F scalar one-row tasks'
    assert accepted['native_first_packets_handled']==accepted['pure_preflights_handled']==3
    assert all(name.startswith('sirius_local_value_') for name in optional['shared_core_functions_observed'])
    action=doc(current/'actions.json')['control-8']
    folder=current/'control-8'
    owner=doc(folder/'owner.json')
    assert owner['passed'] and owner['status']=='terminal' and owner['child_exit']==0
    assert owner['stop_reason'] is None and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['source']==dict(revision=action['revision'],status='')
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152)
    assert action['limits']==owner['limits']==limits
    assert owner['wall_seconds']<=300 and owner['peak_sampled_rss_kib']<=4194304
    assert owner['minimum_sampled_available_kib']>=2097152
    assert accepted['owner']==pin(folder/'owner.json')
    assert accepted['native_observation']==pin(folder/'native-observation.json')
    assert (folder/'inputs-before.json').read_bytes()==(folder/'inputs-after.json').read_bytes()
    assert pin(folder/'owner.json')==manifest['owners'][str((folder/'owner.json').relative_to(root))]['identity']
    build=doc(provider/'provider-build-result.json')
    contract=doc(provider/'source-contract.json')
    proof=doc(provider/'linked-native-proof.json')
    assert build['passed'] and proof['passed']
    assert build['native_objects']==contract['native_objects']==[row['object'] for row in proof['objects']]
    assert len(proof['objects'])==8 and len(proof['raw_v4_objects'])==3
    assert all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
    assert proof['library_identity']==build['library_identity']=={key:manifest['library'][key] for key in ('bytes','sha256')}
    assert proof['library']==build['library']==manifest['library']['path']
    assert proof['component_interface']==accepted['transport_component_scope']['component_interface']
    assert proof['shared_value_interface']['passed'] and proof['shared_value_interface']['core_count']==9
    assert proof['shared_value_interface']['raw_core_calls']==dict(transport=2696,endpoint=3954,dopri_phase=1576)
    for key in ('fixed_signatures_preserved','no_ordinary_writes','no_FP_observers','strict_scope_source_bound'):
        assert proof['shared_value_interface'][key] is True,key
    assert proof['private_value_math_objects']==[proof['shared_value_math_object']]
    assert proof['shared_value_math_object']['all_retained_function_opcodes_match']
    assert sum(map(len,proof['raw_v4_direct_helper_joins'].values()))==8226
    assert proof['native_FMA_mnemonics']==0 and proof['EVEX_confined_to_separately_guarded_v4_functions']
    assert proof['all3_complete_SPIRV_spans_present_once'] and proof['linked_runtime_symbols_resolve']
    assert proof['unchanged_native_indices']==list(range(7))
    assert proof['original_v3_CPU_guard_and_ordinary_fallback_retained']
    for stage in ('transport','endpoint','dopri_phase'):
        table=proof['raw_unwind_tables'][stage]
        for key in ('full_function_FDE','asynchronous','complete_GP_save_rules','shrinkwrap_and_epilogues_checked','exact_relative_CFI_matches_object'):
            assert table[key] is True,(stage,key)
    for key in ('saved_full_mxcsr','traps_masked','round_nearest','daz_ftz_off','restored_full_mxcsr','actual_callback_code_join'):
        assert proof['fp_scope'][key] is True,key
    for stage,row in zip(('transport','endpoint','dopri_phase'),proof['raw_v4_objects']):
        assert row['all_retained_function_opcodes_match']
        assert {'transport':'sirius_closed_transport_packet','endpoint':'sirius_closed_endpoint_packet','dopri_phase':'sirius_closed_dp_packet'}[stage] in {item['symbol'] for item in row['functions']}
    assert optional['provider_linked_proof']==pin(provider/'linked-native-proof.json')
    assert accepted['transport_component_scope']['current_AX_link_proof']==pin(provider/'linked-native-proof.json')
    for key,action_key in (('library','runtime_provider'),('ICD','runtime_icd')):
        path=root/build[key]
        assert str(path)==action[action_key] and pin(path)==build[key+'_identity']
    host=doc(root/'out/native-cpu-owned-lower-increments/host-build-result.json')
    assert host['passed'] and str(root/host['executable'])==action['executable']
    assert pin(Path(action['executable']))==host['executable_identity']
    for name,identity in action['expected_identities'].items():
        assert pin(root/name)==identity,name
    return [root/name for name in manifest['owners']] + [manifest_path,current/'actions.json',current/'control-8-acceptance.json',
        folder/'owner.json',folder/'gtest.xml',folder/'native-observation.json',
        folder/'inputs-before.json',folder/'inputs-after.json',
        provider/'source-contract.json',provider/'provider-build-result.json',
        provider/'linked-native-proof.json',root/'out/native-cpu-owned-lower-increments/host-build-result.json']
