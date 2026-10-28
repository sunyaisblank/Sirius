"""Pinned AU static proof, now exact AX closed3/local9 plus all3 strong FDE joins."""
from pathlib import Path
import hashlib
R=Path.cwd();AUTHORITY=R/'out/native-cpu-shared-value-provider/verify_provider.py'
assert hashlib.sha256(AUTHORITY.read_bytes()).hexdigest()=='55bf38021c9b70292338e9470e4105d7910114c1f8b24fac5e0442eef5904c09'


def replace(text,before,after):
    assert text.count(before)==1,before
    return text.replace(before,after,1)


def prepared_source():
    body=AUTHORITY.read_text().replace('native-cpu-shared-value-provider','native-cpu-closed-all-provider')
    body=body.replace('sirius_shared_value_','sirius_local_value_')
    body=replace(body,"assert contract['shared_value_interface'] == raw_component['shared_value_interface']",
        "assert contract['shared_value_interface'] == raw_component['shared_value_interface']\n"+
        "assert contract['shared_value_interface']['raw_core_calls'] == dict(transport=2696,endpoint=3954,dopri_phase=1576)\n"+
        "assert contract['native_calls'] == raw_component['native_calls']\n"+
        "assert contract['native_receipts'] == raw_component['native_receipts']\n"+
        "assert contract['transport_call_proof'] == raw_component['transport_call_proof']\n"+
        "assert contract['dp_unwind_repair'] == raw_component['dp_unwind_repair']\n"+
        "assert contract['component_baseline_compilations'] == 1")
    body=replace(body,"assert contract['basis_source_identity'] == old['full_candidate_source_identity']",
        "assert contract['basis_source_identity'] == doc(R/'out/native-cpu-closed-dp-unwind-provider/provider/source-contract.json')['full_candidate_source_identity']")
    body=replace(body,"'Sirius-native-cpu-shared-value-v45'","'Sirius-native-cpu-closed-all-v48'")
    body=replace(body,"'SIR-CPU-25.2.8AU'","'SIR-CPU-25.2.8AX'")
    before="""assert len(by_name) == 11 and set(by_name) == {
    'adapter-only.o', 'transport.o', 'endpoint.o', 'dopri_phase.o', 'optional.o', 'dispatch.o', 'v3.o',
    'shared-dp.o', 'shared-transport.o', 'shared-endpoint.cleanup.o', 'shared-math.o'}"""
    after="""assert len(by_name) == 11 and set(by_name) == {
    'adapter-only.o','transport.o','endpoint.o','dopri_phase.o','optional.o','dispatch.o','v3.o',
    'transport.closed.o','endpoint.closed.o','dp.unwind.o','math.o'}"""
    body=replace(body,before,after)
    before="""for item, scope in zip(raw_v4_objects, ('transport', 'endpoint', 'dp')):
    assert any(row['symbol'] == 'sirius_scope_' + scope + '_packet' for row in item['functions'])"""
    after="""for item, scope in zip(raw_v4_objects, ('transport','endpoint','dp')):
    assert len(item['functions']) == 1 and item['functions'][0]['symbol'] == 'sirius_closed_' + scope + '_packet'
assert not [row['name'] for row in linked.symbols if row['section'] and row['name'].startswith(('sirius_shared_value_','sirius_scope_'))]
"""
    body=replace(body,before,after)
    before="'sirius_scope_' + ('dp' if slug == 'dopri_phase' else slug) + '_packet'"
    assert body.count(before)==2;body=body.replace(before,"'sirius_closed_' + ('dp' if slug == 'dopri_phase' else slug) + '_packet'")
    body=replace(body,"n.startswith(('sirius_simd_', 'sirius_scope_'))","n.startswith(('sirius_simd_', 'sirius_closed_'))")
    before="""for item, stage, count in zip(raw_v4_objects, ('transport', 'endpoint', 'dopri_phase'), (2701, 3959, 1581)):
    joins = direct_joins(R / item['object']['path'], item['functions'], lambda n: n.startswith('sirius_local_value_'))
    expected = CORES - ({'sirius_local_value_tagged_sqrt8'} if stage == 'dopri_phase' else set())
    assert len(joins) == count and {row['target'] for row in joins} == expected
    raw_direct_helper_joins[item['object']['path']] = joins
assert sum(map(len, raw_direct_helper_joins.values())) == 8241"""
    after="""from collections import Counter
for item, stage, count in zip(raw_v4_objects, ('transport','endpoint','dopri_phase'), (2696,3954,1576)):
    joins = direct_joins(R/item['object']['path'],item['functions'],lambda n:n.startswith('sirius_local_value_'))
    assert len(joins) == count and dict(Counter(row['target'] for row in joins)) == raw_component['native_calls'][stage]
    raw_direct_helper_joins[item['object']['path']] = joins
assert sum(map(len,raw_direct_helper_joins.values())) == 8226"""
    body=replace(body,before,after)
    anchor='# Count the exact retained functions in these input objects, rather than\n'
    extra="""# Actual object-to-final full FDE/CFI plus GP push/restore machine joins.
# T/E use their own emitted prologue offsets/stack sizes. DP retains the exact
# previous full-function verifier and all of its GP/call/return obligations.
import importlib.util
unwind_path = W/'unwind.py'
checked(unwind_path,contract['source_inputs'][str(unwind_path.relative_to(R))])
spec_unwind = importlib.util.spec_from_file_location('AX_actual_entry_unwind',unwind_path)
unwind_module = importlib.util.module_from_spec(spec_unwind);spec_unwind.loader.exec_module(unwind_module)
unwind_module.Elf = Elf
raw_unwind_tables = {}
for stage,item,symbol in zip(('transport','endpoint','dopri_phase'),raw_v4_objects,
        ('sirius_closed_transport_packet','sirius_closed_endpoint_packet','sirius_closed_dp_packet')):
    source_table = (raw_component['dp_unwind_repair']['object_unwind_table'] if stage=='dopri_phase'
                    else raw_component['native_receipts'][stage]['unwind_table'])
    object_table = unwind_module.verify(R/item['object']['path'],symbol,source_table)
    final_table = unwind_module.verify(library,symbol,source_table)
    for key in ('relative_range','relative_CFI_rows','initial_CFA','body_CFA','saved_GP_rules','saved_return_PC',
                'covered_call_sites','checked_return_offsets'):
        assert final_table[key] == object_table[key],(stage,key)
    for key in ('full_function_FDE','complete_GP_save_rules','shrinkwrap_and_epilogues_checked'):
        assert final_table[key] is True and object_table[key] is True,(stage,key)
    assert final_table['function_bytes'] == item['functions'][0]['bytes']
    if stage=='dopri_phase':
        assert source_table['asynchronous'] is True
    else:
        assert final_table['actual_GP_push_and_restore_machine_bytes_checked'] is True
    final_table.update(asynchronous=True,exact_relative_CFI_matches_object=True,
        exact_object=item['object'],object_unwind_table=object_table,
        verifier=dict(path=str(unwind_path.relative_to(R)),**pin(unwind_path)))
    raw_unwind_tables[stage] = final_table
closed_dp_unwind_table = dict(raw_unwind_tables['dopri_phase'],exact_relative_CFI_matches_repaired_object=True)
"""
    body=replace(body,anchor,extra+anchor)
    body=replace(body,'shared_value_math_object=shared_value_math_object, shared_value_interface=shared_value_interface,',
        "shared_value_math_object=shared_value_math_object, shared_value_interface=shared_value_interface,\n"+
        "    private_value_math_objects=[shared_value_math_object], native_calls=raw_component['native_calls'],\n"+
        "    raw_unwind_tables=raw_unwind_tables, closed_dp_unwind_table=closed_dp_unwind_table,\n"+
        "    dp_unwind_repair=raw_component['dp_unwind_repair'], transport_call_proof=raw_component['transport_call_proof'],")
    start=body.index("    scope='Static AU eleven-input composition only.");end=body.index("')\nwith (P / 'linked-native-proof.json')",start)+2
    body=body[:start]+"    scope='Static AX eleven-part closed3/local9 composition only. First7, both headers/cs_exec_fn, public AS/GCC/dispatcher and baseline/v3/AVX2/scalar paths remain exact. Only private admitted v4 target identifiers select closed canonical T/E; existing closed DP with verified unwind metadata and local-effects cc92 math9 are reused. All8226 actual native PLT32/-4 helper fields bind nine unique hidden core definitions; T2696 accounts for independently audited16-to1 static common-tail joining, E3954 and DP1576 bind exact actual receipts. Every closed object/function/relocation and full function-relative FDE/CFI/GP-save/call/return join is checked at the final DSO. Mixed masks, metadata/phase/root/refusal/transaction/source order and genuine positive-job/full saved-restored FPscope remain required. This is no numerical/performance/scientific/release transfer; fresh original C8/all32 and every original downstream gate remain mandatory.')"+body[end:]
    return body


if __name__=='__main__':
    exec(compile(prepared_source(),str(Path(__file__)),'exec'))
