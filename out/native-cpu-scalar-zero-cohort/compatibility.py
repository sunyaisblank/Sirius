"""Require actual changed scalar Transport execution and exact provider code joins."""
from pathlib import Path
import hashlib
import json


def check(root):
    def doc(path):return json.loads(path.read_text())
    def pin(path):
        with path.open('rb') as stream:
            return {'bytes':path.stat().st_size,'sha256':hashlib.file_digest(stream,'sha256').hexdigest()}
    def named(row):
        path=root/row['path'];assert pin(path)=={k:row[k] for k in ('bytes','sha256')};return path
    work=root/'out/native-cpu-scalar-zero-compatibility'
    for phase in ('build','run'):
        owner=doc(work/phase/'owner.json')
        assert owner['passed'] and owner['status']=='terminal' and owner['child_exit']==0
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    result=doc(work/'result.json');comparison=result['comparison']
    assert result['passed'] and result['runtime_exit']==0 and result['all_prepared_inputs_unchanged']
    for key,value in (('case_count',38),('active_rows',38),('capacity_rows',40),
                      ('packets_per_implementation',5),('raw_calls',10),('compared_output_words',187960)):
        assert comparison[key]==value
    assert comparison['full_word_match'] and comparison['input_unchanged'] and comparison['inactive_sentinels_preserved']
    assert comparison['mismatch_count']==0 and comparison['first_mismatch'] is None
    assert result['output_hashes']['baseline']==result['output_hashes']['candidate']
    provider=root/'out/native-cpu-scalar-zero/provider'
    objects=doc(provider/'source-contract.json')['native_objects'];assert len(objects)==8
    assert result['candidate_transport_object']==objects[0]
    assert result['transport_component_objects']==[objects[0],objects[3]]
    build=doc(named(result['compile_result']));link=doc(named(result['link_proof']))
    assert build['passed'] and build['diagnostic_shim_compilations']==1 and build['new_links']==1
    assert link['passed'] and len(link['objects'])==4
    assert [row['object'] for row in link['objects']]==[result['renamed_baseline_object'],objects[0],objects[3],build['shim_object']]
    assert all(row['all_retained_function_opcodes_match'] for row in link['objects'])
    current=doc(provider/'linked-native-proof.json')
    assert current['passed'] and [row['object'] for row in current['objects']]==objects
    assert all(row['all_retained_function_opcodes_match'] for row in current['objects'])
    for row in objects+[result['renamed_baseline_object'],build['driver_object'],build['shim_object']]:named(row)
    inputs={str(path.relative_to(root)):pin(path) for path in work.rglob('*') if path.is_file()}
    for path in (provider/'source-contract.json',provider/'linked-native-proof.json',named(result['executable'])):
        inputs[str(path.relative_to(root))]=pin(path)
    return {'inputs':inputs,'scope':{'component_result':pin(work/'result.json'),
        'component_link_proof':result['link_proof'],'component_build':result['compile_result'],
        'current_U_link_proof':pin(provider/'linked-native-proof.json'),'compared_output_words':187960,
        'candidate_transport_object':objects[0],
        'scope':'Fresh actual changed scalar-zero Transport, unchanged full38 driver/renamedQbaseline/all187960words. Exact standalone and all8 current-provider opcode joins. One new shim compile/link, zero fixture or baseline math rebuilds. Finite compatibility, no benefit/frame/scientific/release claim.'}}
