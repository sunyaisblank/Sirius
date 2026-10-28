"""Join both sealed finite components to the exact current AT provider."""
from pathlib import Path
import hashlib
import importlib.util
import json


def check(root):
    def pin(path):
        with path.open('rb') as stream:
            return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())
    def doc(path):return json.loads(path.read_text())
    gate_path=root/'out/native-cpu-dispatch-token-provider/component_gate.py'
    spec=importlib.util.spec_from_file_location('AT_sealed_component_gate',gate_path)
    gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
    component=gate.component()
    provider=root/'out/native-cpu-dispatch-token-provider/provider'
    contract=doc(provider/'source-contract.json');proof=doc(provider/'linked-native-proof.json')
    assert proof['passed'] and contract['adapter_input_objects']==component['parts']
    objects=contract['native_objects'];assert len(objects)==8
    assert [row['object'] for row in proof['objects']]==objects
    assert all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
    assert [row['object'] for row in proof['raw_v4_objects']]==component['raw_objects']
    assert all(row['all_retained_function_opcodes_match'] and row['functions'] for row in proof['raw_v4_objects'])
    manifest=doc(component['manifest']);math=doc(component['math_manifest'])
    assert manifest['finite_component_passed'] and not manifest['provider_adopted']
    assert math['finite_component_passed'] and not math['provider_adopted']
    inputs={str(path.relative_to(root)):pin(path) for path in
        (gate_path,component['manifest'],component['archive'],component['math_manifest'],component['math_archive'],
         provider/'source-contract.json',provider/'linked-native-proof.json')}
    for row in component['parts']+objects:
        path=root/row['path'];identity={key:row[key] for key in ('bytes','sha256')}
        assert pin(path)==identity;inputs[row['path']]=identity
    return dict(inputs=inputs,scope=dict(component_result=pin(component['manifest']),
        component_link_proof=manifest['members']['out/native-cpu-dispatch-admission-token/transport-link-proof.json'],
        component_compile_result=manifest['members']['out/native-cpu-dispatch-admission-token/build-result.json'],
        math_component_result=pin(component['math_manifest']),current_AT_link_proof=pin(provider/'linked-native-proof.json'),
        candidate_transport_object=component['raw_objects'][0],compared_output_words=187960,
        scope='Sealed original R38/cap40/all187960-word Transport and focused Endpoint fixtures each ran once in public and admitted modes with every original assertion/count retained. Private diagnostic calls follow actual unchanged pure preflight; failed preflight invokes no trusted worker. Current AT binds new baseline scalar3/adapter, unchanged pure admission, ASraw3/AVX2 fallback and math3. Current C8 must independently observe the actual positive live job and hidden workers; no predecessor acceptance or benefit/frame/scientific/release transfer.'))
