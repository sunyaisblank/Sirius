"""Require fresh AU C8 and exact current static composition before any ABBA slot."""
from pathlib import Path
import hashlib
import json


def check():
    root=Path.cwd()
    current=root/'out/native-cpu-shared-value-cohort'
    provider=root/'out/native-cpu-shared-value-provider/provider'
    def doc(path):
        return json.loads(path.read_text())
    def pin(path):
        with path.open('rb') as stream:
            return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())
    accepted=doc(current/'control-8-acceptance.json')
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
    build=doc(provider/'provider-build-result.json')
    contract=doc(provider/'source-contract.json')
    proof=doc(provider/'linked-native-proof.json')
    assert build['passed'] and proof['passed']
    assert build['native_objects']==contract['native_objects']==[row['object'] for row in proof['objects']]
    assert len(proof['objects'])==8 and len(proof['raw_v4_objects'])==3
    assert all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
    for stage,row in zip(('transport','endpoint','dopri_phase'),proof['raw_v4_objects']):
        assert row['all_retained_function_opcodes_match']
        assert {'transport':'sirius_scope_transport_packet','endpoint':'sirius_scope_endpoint_packet','dopri_phase':'sirius_scope_dp_packet'}[stage] in {item['symbol'] for item in row['functions']}
    assert optional['provider_linked_proof']==pin(provider/'linked-native-proof.json')
    assert accepted['transport_component_scope']['current_AU_link_proof']==pin(provider/'linked-native-proof.json')
    for key,action_key in (('library','runtime_provider'),('ICD','runtime_icd')):
        path=root/build[key]
        assert str(path)==action[action_key] and pin(path)==build[key+'_identity']
    host=doc(root/'out/native-cpu-owned-lower-increments/host-build-result.json')
    assert host['passed'] and str(root/host['executable'])==action['executable']
    assert pin(Path(action['executable']))==host['executable_identity']
    for name,identity in action['expected_identities'].items():
        assert pin(root/name)==identity,name
    return [current/'actions.json',current/'control-8-acceptance.json',
        folder/'owner.json',folder/'gtest.xml',folder/'native-observation.json',
        folder/'inputs-before.json',folder/'inputs-after.json',
        provider/'source-contract.json',provider/'provider-build-result.json',
        provider/'linked-native-proof.json',root/'out/native-cpu-owned-lower-increments/host-build-result.json']
