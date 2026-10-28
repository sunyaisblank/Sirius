"""Bind the mixed AV static interface; C8 supplies its first real numerical gate."""
from pathlib import Path
import hashlib,importlib.util,json


def check(root):
 def pin(p):
  with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())
 def doc(p):return json.loads(p.read_text())
 gate_path=root/'out/native-cpu-closed-dp-unwind-provider/component_gate.py'
 spec=importlib.util.spec_from_file_location('AV_compiled_gate',gate_path);gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
 component=gate.component()
 expected=dict(compiled_interface_passed=True,static_symbol_binding_passed=True,numerical_qualification_transferred=False,genuine_job_numerical_gate_required=True)
 assert component['component_interface']==expected
 provider=root/'out/native-cpu-closed-dp-unwind-provider/provider';contract=doc(provider/'source-contract.json');proof=doc(provider/'linked-native-proof.json')
 assert proof['passed'] and contract['adapter_input_objects']==component['parts']
 assert contract['dp_unwind_repair']==proof['dp_unwind_repair']==component['dp_unwind_repair']
 table=proof['closed_dp_unwind_table'];object_table=component['dp_unwind_repair']['object_unwind_table']
 for key in ('full_function_FDE','asynchronous','complete_GP_save_rules','shrinkwrap_and_epilogues_checked','exact_relative_CFI_matches_repaired_object'):assert table[key] is True,key
 assert table['relative_range']==object_table['relative_range'] and table['relative_CFI_rows']==object_table['relative_CFI_rows']
 objects=contract['native_objects'];assert len(objects)==8 and len(component['parts'])==12
 assert [row['object'] for row in proof['objects']]==objects
 assert all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
 assert [row['object'] for row in proof['raw_v4_objects']]==component['raw_v4_objects']
 names=('sirius_scope_transport_packet','sirius_scope_endpoint_packet','sirius_closed_dp_packet')
 for name,row in zip(names,proof['raw_v4_objects']):
  assert row['all_retained_function_opcodes_match'] and name in {item['symbol'] for item in row['functions']}
 suffixes=('dense_add8','dense_multiply8','dense_divide8','dense_sqrt8','tagged_add8','tagged_subtract8','tagged_multiply8','tagged_divide8','tagged_sqrt8')
 for key,prefix in [('shared_value_math_object','sirius_shared_value_'),('closed_dp_math_object','sirius_local_value_')]:
  item=proof[key];assert item['object']==component[key] and item['all_retained_function_opcodes_match']
  assert {prefix+s for s in suffixes}<={row['symbol'] for row in item['functions']}
 interface=proof['shared_value_interface'];closed=proof['closed_dp_interface']
 assert interface['passed'] and interface['core_count']==9 and interface['raw_core_calls']==dict(transport=2701,endpoint=3959)
 assert closed['passed'] and closed['core_count']==9 and closed['raw_core_calls']==dict(dopri_phase=1576)
 for item in [interface,closed]:
  for key in ['fixed_signatures_preserved','no_ordinary_writes','no_FP_observers','strict_scope_source_bound']:assert item[key] is True,key
 for key in ['saved_full_mxcsr','traps_masked','round_nearest','daz_ftz_off','restored_full_mxcsr','actual_callback_code_join']:assert proof['fp_scope'][key] is True,key
 assert proof['callback_live_job_register']['no_intervening_job_register_writes']
 assert proof['actual_hidden_admitted_worker_calls'] and proof['job_observation_layout']
 compiled=root/'out/native-cpu-closed-dp-unwind-component/build-result.json';assert doc(compiled)['passed']
 inputs={str(p.relative_to(root)):pin(p) for p in [gate_path,compiled,provider/'source-contract.json',provider/'linked-native-proof.json']}
 for row in component['parts']+objects:
  p=root/row['path'];identity={k:row[k] for k in ['bytes','sha256']};assert pin(p)==identity;inputs[row['path']]=identity
 return dict(inputs=inputs,scope=dict(component_result=pin(compiled),component_compile_result=pin(compiled),component_interface=expected,
  current_AW_link_proof=pin(provider/'linked-native-proof.json'),candidate_transport_object=component['raw_v4_objects'][0],
  scope='AW current original numerical C8/all32 for same positive-job closed DP/local9 machine code with repaired unwind metadata; unchanged AU T/E/shared9 and all prior public/scalar/v3/AVX2 paths remain independently bound. Current static checks do not transfer numerical/benefit/frame/scientific/release acceptance.'))
