"""Bind the compiled private interface; numerical acceptance starts at AU C8."""
from pathlib import Path
import hashlib
import importlib.util
import json


def check(root):
    def pin(path):
        with path.open('rb') as stream:
            return dict(bytes=path.stat().st_size,
                        sha256=hashlib.file_digest(stream,'sha256').hexdigest())
    def doc(path):
        return json.loads(path.read_text())
    gate_path=root/'out/native-cpu-shared-value-provider/component_gate.py'
    spec=importlib.util.spec_from_file_location('AU_compiled_interface_gate',gate_path)
    gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
    component=gate.component()
    expected_interface=dict(compiled_interface_passed=True,static_symbol_binding_passed=True,
        numerical_qualification_transferred=False,genuine_job_numerical_gate_required=True)
    assert component['component_interface']==expected_interface
    provider=root/'out/native-cpu-shared-value-provider/provider'
    contract=doc(provider/'source-contract.json');proof=doc(provider/'linked-native-proof.json')
    assert proof['passed'] and contract['adapter_input_objects']==component['parts']
    objects=contract['native_objects'];assert len(objects)==8 and len(component['parts'])==11
    assert [row['object'] for row in proof['objects']]==objects
    assert all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
    assert [row['object'] for row in proof['raw_v4_objects']]==component['raw_v4_objects']
    names=('sirius_scope_transport_packet','sirius_scope_endpoint_packet','sirius_scope_dp_packet')
    assert len(proof['raw_v4_objects'])==3
    for name,row in zip(names,proof['raw_v4_objects']):
        assert row['all_retained_function_opcodes_match'] and row['functions']
        assert name in {item['symbol'] for item in row['functions']}
    math=proof['shared_value_math_object']
    assert math['object']==component['shared_value_math_object']
    assert math['all_retained_function_opcodes_match']
    cores={'sirius_shared_value_'+name for name in (
        'dense_add8','dense_multiply8','dense_divide8','dense_sqrt8',
        'tagged_add8','tagged_subtract8','tagged_multiply8','tagged_divide8','tagged_sqrt8')}
    assert cores<={row['symbol'] for row in math['functions']}
    interface=proof['shared_value_interface']
    assert interface['passed'] and interface['core_count']==9
    assert interface['raw_core_calls']==dict(dopri_phase=1581,transport=2701,endpoint=3959)
    for key in ('fixed_signatures_preserved','no_ordinary_writes','no_FP_observers','strict_scope_source_bound'):
        assert interface[key] is True,key
    for key in ('saved_full_mxcsr','traps_masked','round_nearest','daz_ftz_off',
                'restored_full_mxcsr','actual_callback_code_join'):
        assert proof['fp_scope'][key] is True,key
    assert proof['callback_live_job_register']['no_intervening_job_register_writes']
    assert proof['actual_hidden_admitted_worker_calls'] and proof['job_observation_layout']
    compiled=root/'out/native-cpu-shared-value-component/build-result.json'
    assert doc(compiled)['passed']
    inputs={str(path.relative_to(root)):pin(path) for path in
        (gate_path,compiled,provider/'source-contract.json',provider/'linked-native-proof.json')}
    for row in component['parts']+objects:
        path=root/row['path'];identity={key:row[key] for key in ('bytes','sha256')}
        assert pin(path)==identity;inputs[row['path']]=identity
    return dict(inputs=inputs,scope=dict(component_result=pin(compiled),
        component_compile_result=pin(compiled),component_interface=expected_interface,
        current_AU_link_proof=pin(provider/'linked-native-proof.json'),
        candidate_transport_object=component['raw_v4_objects'][0],
        scope='Compiled/static shared-value interface only; no finite numerical pass is asserted or transferred. AU C8/all32 on the unchanged H host must first execute the actual positive-job private scoped entries and shared cores inside the complete saved/masked/restored MXCSR callback. Original AT public wrappers/AS raw3/GCC math and scalar/AVX2/v3 fallback stay separately preserved. No predecessor numerical/benefit/frame/scientific/release acceptance transfers.'))
