"""AX compile/static closure only; the unchanged original C8 is its numerical gate."""
from pathlib import Path
import hashlib,importlib.util,json


def check(root):
 def pin(p):
  with p.open('rb') as s:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(s,'sha256').hexdigest())
 def doc(p):return json.loads(p.read_text())
 gate_path=root/'out/native-cpu-closed-all-provider/component_gate.py'
 spec=importlib.util.spec_from_file_location('AX_compiled_gate',gate_path);gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
 component=gate.component()
 expected=dict(compiled_interface_passed=True,static_symbol_binding_passed=True,numerical_qualification_transferred=False,genuine_job_numerical_gate_required=True)
 assert component['component_interface']==expected
 provider=root/'out/native-cpu-closed-all-provider/provider'
 contract=doc(provider/'source-contract.json');proof=doc(provider/'linked-native-proof.json')
 assert proof['passed'] and contract['adapter_input_objects']==component['parts'] and len(component['parts'])==11
 objects=contract['native_objects'];assert len(objects)==8 and [row['object'] for row in proof['objects']]==objects
 assert all(row['all_retained_function_opcodes_match'] for row in proof['objects'])
 assert [row['object'] for row in proof['raw_v4_objects']]==component['raw_v4_objects']
 for name,row in zip(('sirius_closed_transport_packet','sirius_closed_endpoint_packet','sirius_closed_dp_packet'),proof['raw_v4_objects']):
  assert row['all_retained_function_opcodes_match'] and len(row['functions'])==1 and row['functions'][0]['symbol']==name
 math=proof['shared_value_math_object'];assert math['object']==component['shared_value_math_object'] and math['all_retained_function_opcodes_match']
 suffixes=('dense_add8','dense_multiply8','dense_divide8','dense_sqrt8','tagged_add8','tagged_subtract8','tagged_multiply8','tagged_divide8','tagged_sqrt8')
 assert {'sirius_local_value_'+s for s in suffixes}<={row['symbol'] for row in math['functions']}
 assert proof['private_value_math_objects']==[math]
 interface=proof['shared_value_interface'];assert interface['passed'] and interface['core_count']==9
 assert interface['raw_core_calls']==dict(transport=2696,endpoint=3954,dopri_phase=1576)
 for key in ('fixed_signatures_preserved','no_ordinary_writes','no_FP_observers','strict_scope_source_bound'):assert interface[key] is True,key
 assert proof['native_calls']==component['native_calls'] and sum(map(len,proof['raw_v4_direct_helper_joins'].values()))==8226
 for stage in ('transport','endpoint','dopri_phase'):
  table=proof['raw_unwind_tables'][stage]
  for key in ('full_function_FDE','asynchronous','complete_GP_save_rules','shrinkwrap_and_epilogues_checked','exact_relative_CFI_matches_object'):assert table[key] is True,(stage,key)
  assert table['relative_range']==table['object_unwind_table']['relative_range']
  assert table['relative_CFI_rows']==table['object_unwind_table']['relative_CFI_rows']
  if stage!='dopri_phase':assert table['actual_GP_push_and_restore_machine_bytes_checked'] is True
 dp=proof['closed_dp_unwind_table'];assert dp['exact_relative_CFI_matches_repaired_object'] is True
 assert proof['dp_unwind_repair']==component['dp_unwind_repair']==contract['dp_unwind_repair']
 for key in ('saved_full_mxcsr','traps_masked','round_nearest','daz_ftz_off','restored_full_mxcsr','actual_callback_code_join'):assert proof['fp_scope'][key] is True,key
 assert proof['callback_live_job_register']['no_intervening_job_register_writes']
 assert proof['actual_hidden_admitted_worker_calls'] and proof['job_observation_layout']
 compiled=root/'out/native-cpu-closed-all-component/build-result.json';assert doc(compiled)['passed']
 inputs={str(p.relative_to(root)):pin(p) for p in (gate_path,compiled,provider/'source-contract.json',provider/'linked-native-proof.json')}
 for row in component['parts']+objects:
  p=root/row['path'];identity={k:row[k] for k in ('bytes','sha256')};assert pin(p)==identity;inputs[row['path']]=identity
 return dict(inputs=inputs,scope=dict(component_result=pin(compiled),component_compile_result=pin(compiled),component_interface=expected,
  current_AX_link_proof=pin(provider/'linked-native-proof.json'),candidate_transport_object=component['raw_v4_objects'][0],
  scope='AX current original numerical C8/all32 for actual closed canonical T/E/DP and existing local-effects register math9 inside genuine positive-job/full FPscope. Public/scalar/GCC/v3/AVX2 fallback remains separately exact. All3 final FDE/CFI/GP joins are prerequisites, not runtime unwind acceptance. No prior numerical/benefit/frame/scientific/release transfer.'))
