"""Compile-interface source closure for AX closed T/E/DP; no numerical transfer."""
from pathlib import Path
import hashlib
import importlib.util
import json

R = Path.cwd()
C = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-closed-all-component'
PREVIOUS = R / 'out/native-cpu-closed-dp-unwind-component'
REMAINING = R / 'out/native-cpu-canonical-closed-remaining'
V = R / 'out/native-cpu-strict-ir-linking'
path = PREVIOUS / 'component.py'
# Executed predecessor sources are bound by their owned immutable build plan.
expected = json.loads((PREVIOUS / 'build-plan.json').read_text())['inputs'][str(path.relative_to(R))]
assert path.stat().st_size == expected['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest() == expected['sha256']
spec = importlib.util.spec_from_file_location('AX_exact_AW_component', path)
previous = importlib.util.module_from_spec(spec); spec.loader.exec_module(previous)
doc, pin, checked, terminal, contract_inputs = previous.doc, previous.pin, previous.checked, previous.terminal, previous.contract_inputs
HEADERS, LOCAL_CORES, DP_CALLS = previous.HEADERS, previous.LOCAL_CORES, previous.DP_CALLS
CORES = LOCAL_CORES
RAW_NAMES = ('sirius_closed_transport_packet', 'sirius_closed_endpoint_packet', 'sirius_closed_dp_packet')
ADAPTER_BASIS = R / 'out/native-cpu-closed-dp-component/rows_adapter.cpp'
ADAPTER_BASIS_PIN = dict(bytes=9093, sha256='dcb841d105eea00200bdadad9d818efffb31eff1980ea3e8c0e8f871388fa2b5')
TRANSPORT_OBJECT = dict(path='bin/linux-gcc/native-cpu-canonical-closed-transport/transport.closed.o',
    bytes=3571376, sha256='c98d747b85ed217df3884cfb50c0bd419be4e90845d0b3d9248d2f5a098d1473')
DP_OBJECT = dict(path='bin/linux-gcc/native-cpu-closed-dp-unwind/dp.unwind.o',bytes=2851176,
    sha256='fd46489c0b6603dcc323cd4ca62206492a02c1182f95a7e98315ddb78a394a24')


def adapter_source():
    checked(ADAPTER_BASIS, ADAPTER_BASIS_PIN)
    before = ADAPTER_BASIS.read_text(); after = before
    for stage in ('transport','endpoint'):
        a,b = 'sirius_scope_'+stage+'_packet', 'sirius_closed_'+stage+'_packet'
        assert after.count(a)==2 and b not in after
        after=after.replace(a,b)
    restored=after
    for stage in ('transport','endpoint'):
        restored=restored.replace('sirius_closed_'+stage+'_packet','sirius_scope_'+stage+'_packet')
    assert restored==before
    prefix=before.index('extern "C" __attribute__((visibility("hidden"))) void')
    assert after[:prefix]==before[:prefix]
    assert (C/'rows_adapter.cpp').read_text()==after
    derivation=doc(C/'derivation.json')
    assert derivation['source_inverse_exact'] and derivation['public_wrapper_prefix_byte_preserved']
    checked(C/'rows_adapter.cpp',derivation['derived'])
    assert derivation['source_inputs'][str(ADAPTER_BASIS.relative_to(R))]==ADAPTER_BASIS_PIN
    return derivation


def native(stage):
    result=doc(REMAINING/(stage+'-native-result.json'))
    assert result['phase']=='native' and result['stage']==stage and result['numerical_qualification_transferred'] is False
    assert result['sole_hidden_C_ABI_entry']=='sirius_closed_'+stage+'_packet'
    assert not result['emitted_fused_instructions'] and result['instructions']>0
    output=result['output']; assert output['path']=='bin/linux-gcc/native-cpu-canonical-closed-'+stage+'/'+stage+'.closed.o'
    checked(R/output['path'],output)
    for phase in ('context-plugin','context','cleanup','native'):
        receipt=doc(REMAINING/(stage+'-'+phase+'-result.json'))
        owner=receipt.get('accepted_artifact_owner',str((REMAINING/(phase+'-'+stage)/'owner.json').relative_to(R)))
        if stage=='endpoint' and phase in ('cleanup','native'):
            expected={'cleanup':'cleanup-canonical-audit-endpoint','native':'native-audited-endpoint'}[phase]
            assert owner==str((REMAINING/expected/'owner.json').relative_to(R))
        terminal(R/owner)
    cleanup=doc(REMAINING/(stage+'-cleanup-result.json'))
    assert cleanup['owned_context_allocas']==0 and cleanup['definitions']==1
    expected_words={'transport':37592,'endpoint':30632}[stage]
    assert cleanup['original_transactional_staging_bytes']==expected_words*4
    staged='%staged_output = alloca ['+str(expected_words)+' x i32], align 64'
    assert cleanup['allocas'].count(staged)==1
    other=[row for row in cleanup['allocas'] if row!=staged]
    if stage=='endpoint':
        audit_path=R/cleanup['residual_audit']['path'];checked(audit_path,cleanup['residual_audit']);audit=doc(audit_path)
        for key in ('passed','read_only','full_module_fingerprint_unchanged','FP_call_inventory_unchanged',
                    'text_reference_predecessor_comments_normalized','strict_dynamic_exceptions',
                    'ordinary_entry_async_uwtable','math9_fixed_cc92_inaccessibleRW'):assert audit[key] is True,key
        assert cleanup['cleanup_repeated'] is False and cleanup['original_failed_owner']==pin(REMAINING/'cleanup-endpoint/owner.json')
        assert audit['residual_bytes']==320 and audit['residual_memory_operations']==44 and audit['value_call_total']==3954
        assert audit['staging_bytes']==122528 and audit['owned_context_allocas']==0
        assert len(other)==6 and {row['name'] for row in audit['residuals']}=={
            'programBase_3.i','count_2.i','layers_2.i','prefix_5.i','component_4.i','root_1.i'}
        assert all(row['full_initialization_dominates_every_memory_use'] and row['no_capture_or_dynamic_address'] for row in audit['residuals'])
    else:
        assert len(other)<=1 and all(('alloca [125 x <8 x i32>], align 64' in row or 'alloca [4000 x i8], align 64' in row) for row in other)
    assert result['cleaned_IR_value_calls']==cleanup['value_calls']
    table=result['unwind_table']; assert table['full_function_FDE'] and table['relative_CFI_rows']
    assert table['relative_range']==[0,table['function_bytes']] and table['covered_call_sites']>0 and table['checked_return_offsets']
    assert set(result['core_relocations'])<=CORES and all(n>0 for n in result['core_relocations'].values())
    if stage=='endpoint':
        assert output==dict(path='bin/linux-gcc/native-cpu-canonical-closed-endpoint/endpoint.closed.o',
            bytes=4846432,sha256='52820e2af329dcce5bb4f84e86138b305db1f92828b05ba63eb4ad9bf558b427')
        assert sum(result['core_relocations'].values())==3954
        assert result['core_relocations']==result['cleaned_IR_value_calls']
    if stage=='transport':
        assert output==TRANSPORT_OBJECT and sum(result['core_relocations'].values())==2696
        assert result['core_relocations']['sirius_local_value_dense_divide8']==1
        assert result['cleaned_IR_value_calls']['sirius_local_value_dense_divide8']==16
    return result


def ingredients():
    old=previous.component(); derivation=adapter_source()
    transport,endpoint=native('transport'),native('endpoint')
    assert old['raw_v4_objects'][2]==DP_OBJECT
    checked(R/DP_OBJECT['path'],DP_OBJECT)
    assert old['closed_dp_math_object']['path']=='bin/linux-gcc/native-cpu-local-effects/math.o'
    products=old['parts'][1:7]+[transport['output'],endpoint['output'],DP_OBJECT,old['closed_dp_math_object']]
    assert len(products)==10
    for row in products: checked(R/row['path'],row)
    inputs=dict(old['source_inputs']); inputs.update(derivation['source_inputs'])
    plans=[]
    for stage in ('transport','endpoint'):
        plan_path=REMAINING/(stage+'-context-plan.json'); plan=doc(plan_path)
        inputs.update(plan['inputs']); plans.append(plan_path)
        for phase in ('context-plugin','context','cleanup','native'):
            for name in (stage+'-'+phase+'-result.json',phase+'-'+stage+'/owner.json',
                         phase+'-'+stage+'/inputs-before.json',phase+'-'+stage+'/inputs-after.json'):
                p=REMAINING/name
                if name.endswith('/owner.json') or '/inputs-' in name:
                    phase=name.split('-'+stage+'/')[0]
                    receipt=doc(REMAINING/(stage+'-'+phase+'-result.json'))
                    accepted=receipt.get('accepted_artifact_owner')
                    if accepted:
                        parent=(R/accepted).parent;p=parent/p.name
                inputs[str(p.relative_to(R))]=pin(p)
    audit_plan=REMAINING/'endpoint-canonical-residual-plan.json'
    inputs.update(doc(audit_plan)['inputs']);inputs[str(audit_plan.relative_to(R))]=pin(audit_plan)
    accepted_cleanup=doc(REMAINING/'endpoint-cleanup-result.json')
    for p in (R/accepted_cleanup['residual_audit']['path'],REMAINING/'endpoint-native-audited-plan.json'):
        inputs[str(p.relative_to(R))]=pin(p)
    inputs.update(doc(REMAINING/'endpoint-native-audited-plan.json')['inputs'])
    # Actual backend-tail merge evidence is an additional adoption prerequisite;
    # its fixed reader is supplied by root after the independent CFG audit.
    audit_path=REMAINING/'transport-native-call-proof.json'
    audit=doc(audit_path)
    assert audit['passed'] and audit['numerical_qualification_transferred'] is False
    assert audit['object']==transport['output']
    assert audit['IR']==doc(REMAINING/'transport-cleanup-result.json')['ir']
    checked(R/audit['IR']['path'],audit['IR'])
    assert audit['IR_call_count']==2711 and audit['native_relocation_count']==2696
    assert audit['dense_divide_IR_sites']==16 and audit['dense_divide_native_targets']==1
    assert audit['all_16_original_call_paths_reach_the_same_exact_target'] is True
    inputs[str(audit_path.relative_to(R))]=pin(audit_path)
    for p in [path,ADAPTER_BASIS,PREVIOUS/'build-plan.json',PREVIOUS/'build-result.json',*plans,*C.glob('*.py'),C/'rows_adapter.cpp',C/'derivation.json']:
        inputs[str(p.relative_to(R))]=pin(p)
    for row in products: inputs[row['path']]={k:row[k] for k in ('bytes','sha256')}
    for name,expected in inputs.items(): checked(R/name,expected)
    calls={'transport':transport['core_relocations'],'endpoint':endpoint['core_relocations'],'dopri_phase':DP_CALLS}
    interface=dict(passed=True,core_count=9,raw_core_calls={stage:sum(rows.values()) for stage,rows in calls.items()},
        fixed_signatures_preserved=True,no_ordinary_writes=True,no_FP_observers=True,strict_scope_source_bound=True)
    assert interface['raw_core_calls']==dict(transport=2696,endpoint=3954,dopri_phase=1576)
    return dict(parent=old,products=products,inputs=inputs,derivation=derivation,
        raw_v4_objects=[transport['output'],endpoint['output'],DP_OBJECT],native_calls=calls,
        native_receipts=dict(transport=transport,endpoint=endpoint),shared_value_math_object=old['closed_dp_math_object'],
        shared_value_interface=interface,dp_unwind_repair=old['dp_unwind_repair'],transport_call_proof=audit)


def component():
    value=ingredients(); terminal(C/'build/owner.json')
    plan,result=doc(C/'build-plan.json'),doc(C/'build-result.json')
    parts=[result['adapter']]+value['products']
    assert len(parts)==11 and plan['products']==value['products'] and result['parts']==parts
    assert result['passed'] and result['compiled_interface_passed'] and result['all_raw_scope_and_core_references_bound']
    assert result['new_native_compilations']==plan['new_native_compilations']==1
    assert result['numerical_qualification_transferred'] is False
    assert set(result['core_definitions'])==CORES and len(result['core_definitions'])==9
    assert result['closed_raw_exports']==list(RAW_NAMES)
    assert plan['merge_argv']==['/usr/bin/ld','-r','-o',str(B/'rows_adapter.o'),str(B/'adapter-only.o')]+[str(R/row['path']) for row in value['products']]
    checked(R/result['adapter']['path'],result['adapter']);checked(R/result['output']['path'],result['output'])
    assert result['output']['path']==str((B/'rows_adapter.o').relative_to(R))
    inputs=dict(value['inputs'])
    for name in ('build-plan.json','build-result.json','build/owner.json','build/inputs-before.json','build/inputs-after.json'):
        p=C/name;inputs[str(p.relative_to(R))]=pin(p)
    old=value['parent']
    return dict(old,merged=result['output'],parts=parts,source_inputs=inputs,
        raw_v4_objects=value['raw_v4_objects'],shared_value_math_object=value['shared_value_math_object'],
        shared_value_interface=value['shared_value_interface'],native_calls=value['native_calls'],
        native_receipts=value['native_receipts'],transport_call_proof=value['transport_call_proof'],
        dp_unwind_repair=value['dp_unwind_repair'],
        component_interface=dict(compiled_interface_passed=True,static_symbol_binding_passed=True,
            numerical_qualification_transferred=False,genuine_job_numerical_gate_required=True),
        scope='AX static11-part compiled interface only: one baseline adapter changes exactly private T/E v4 target identifiers; six old public AS/GCC/dispatcher objects, first7 native rows, DP unwind object and local-effects math9 remain exact. New closed T/E preserve mixed masks, full transactional staging, source/FP/order/contracts. No numerical acceptance is transferred; current genuine-job original C8/all32 is required.')
