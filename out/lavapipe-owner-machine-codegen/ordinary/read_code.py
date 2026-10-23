"""Read saved three finite captures; no subprocesses, dispatches or JIT calls."""
from pathlib import Path
import hashlib, json
R = Path.cwd(); W = Path(__file__).resolve().parent
def doc(p): return json.loads(p.read_text())
def pin(p):
    with p.open('rb') as stream: digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes':p.stat().st_size,'sha256':digest}
def ir_bytes(folder, verification):
    export = verification['ordinary_ir_export']; p = folder/export['path']
    assert pin(p) == {k:export[k] for k in ['bytes','sha256']}
    return p, p.read_bytes()
compiler = doc(R/'out/native-cooperative-transport-intergroup/compiler-result.json')
provider = doc(W.parent/'provider-build-result.json')
assert provider['passed']
driver = R/provider['library']; assert pin(driver) == provider['library_identity']
host = R/'bin/linux-gcc/native-cooperative-transport-intergroup/matched/audited_backend_tests'
assert pin(host)['sha256'] == 'f958a60563b64c79d3f2fb6da2815bde0effdbbe49f121034cfd212ab17b27ff'
rows = {}
for label, kind in [('transport','Transport'),('endpoint','Endpoint'),('dopri_phase','DopriPhase')]:
    folder = W/label; owner = doc(folder/'owner.json')
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['stop_reason'] is None and not owner['cleanup_errors']
    assert owner['remaining_owned_processes'] == []
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
    assert doc(folder/'inputs-before.json') == doc(folder/'inputs-after.json')
    verification = doc(folder/'capture-verification.json'); assert verification['passed']
    assert verification['actual_provider_sha256'] == provider['library_identity']['sha256']
    assert verification['actual_complete_module']['sha256'] == compiler['modules'][kind]['sha256']
    entry = verification['actual_entry']; assert entry['dispatch_id'] == 0
    assert entry['entry_sysv_integer_arguments_2_3_4'] == [64,1,1]
    assert entry['actual_mxcsr'] & 0xe040 == 0 and verification['callback_restoration_observed']
    bindings = doc(folder/'bindings.json')
    assert len(bindings['dispatches']) == 1 and bindings['dispatches'][0]['completed']
    assert bindings['dispatches'][0]['groups_x'] == 1
    event, = [x for x in bindings['jit_events'] if x.get('file') == entry['file']]
    assert str(driver) in event['maps_at_first_call']
    assert pin(folder/event['file'])['sha256'] == event['sha256']
    for section in event['live_sections_after_pointer_lookup']:
        assert pin(folder/section['path']) == {k:section[k] for k in ['bytes','sha256']}
    selected = [x for x in bindings['pointer_lookups'] if entry['file'] in x['matched_files']]
    generation, = {x['module_generation'] for x in selected}
    pointer, = {x['pointer'] for x in selected}
    assert pointer == entry['actual_entry_pc']
    policy = doc(folder/'machine-policy.json')
    lookups = [x for x in policy['lookups'] if x['completed'] and
               x['module_generation'] == generation and x['pointer'] == pointer]
    assert lookups
    engine, = {x['engine'] for x in lookups}
    creation, = [x for x in policy['engines'] if x['completed'] and
                  x['engine'] == engine and x['module_generation'] == generation]
    assert creation['initial_OptLevel'] == 2
    emissions = [x for x in policy['emissions'] if x['completed'] and
                 x['engine'] == engine and x['module_generation'] == generation]
    assert emissions and all(any(l['sequence'] < e['sequence'] < e['return_sequence'] < l['return_sequence']
           and l['thread'] == e['thread'] for l in lookups) for e in emissions)
    emitted_ids = {x['id'] for x in emissions}
    factories = [x for x in policy['factories'] if x['emission_id'] in emitted_ids]
    assert factories and all(x['engine'] == engine and x['module_generation'] == generation for x in factories)
    for factory in factories:
        emission = policy['emissions'][factory['emission_id']]
        assert emission['sequence'] < factory['sequence'] < emission['return_sequence']
        assert factory['thread'] == emission['thread']
    fast = sum(x['allocator'] == 'fast' for x in factories)
    greedy = sum(x['allocator'] == 'greedy' for x in factories)
    bridges = [x for x in policy['bridges'] if x['engine'] == engine and x['module_generation'] == generation]
    if kind == 'Transport':
        assert not bridges and greedy > 0 and fast == 0
    else:
        bridge, = bridges
        assert bridge['completed'] and bridge['returned_bool'] == 1
        assert creation['return_sequence'] < bridge['sequence'] < bridge['return_sequence']
        assert all(bridge['return_sequence'] < e['sequence'] for e in emissions)
        assert fast > 0 and greedy == 0
        completed = verification['ordinary_ir_export']['successful_optimizer_pass_ids']
        assert set(completed) <= set(bridge['completed_optimizer_pass_ids_at_bridge'])
    prior_root = R/('out/native-original-input-ownership/ordinary' if kind == 'Transport'
                   else 'out/lavapipe-intergroup-row-owner/ordinary')
    prior_folder = prior_root/label
    prior_v = doc(prior_folder/'capture-verification.json'); assert prior_v['passed']
    assert prior_v['actual_complete_module']['sha256'] == compiler['modules'][kind]['sha256']
    current_path, current_ir = ir_bytes(folder, verification)
    prior_path, prior_ir = ir_bytes(prior_folder, prior_v)
    # Deliberately no text/SSA/address normalization. A mismatch refuses this
    # exact-IR gate rather than silently broadening the equivalence claim.
    assert current_ir == prior_ir, ('Optimized IR changed', kind, pin(current_path), pin(prior_path))
    prior_review = doc(prior_root/'ordinary-code-review.json')
    assert prior_review['structural_premise_observed']
    prior_row = prior_review['rows'][kind]
    rows[kind] = {'structural_premise_observed':True,
        'machine_policy_isolation_observed':True, 'optimized_IR_unchanged':True,
        'ordinary_code_source_proof_reused':str((prior_root/'ordinary-code-review.json').relative_to(R)),
        'prior_proof_identity':pin(prior_root/'ordinary-code-review.json'),
        'actual_complete_module':verification['actual_complete_module'],
        'actual_capture_verification':pin(folder/'capture-verification.json'),
        'actual_optimized_IR':{'path':str(current_path.relative_to(R)),**pin(current_path)},
        'prior_optimized_IR':{'path':str(prior_path.relative_to(R)),**pin(prior_path)},
        'actual_engine':engine, 'actual_module_generation':generation,
        'initial_OptLevel':creation['initial_OptLevel'],
        'post_IR_bridge_successes':len(bridges),
        'fast_allocator_factories_in_selected_engine_emission':fast,
        'greedy_allocator_factories_in_selected_engine_emission':greedy,
        'actual_emission_envelopes':emissions,
        'actual_allocator_factories':factories,
        'actual_original_block_arguments':entry['entry_sysv_integer_arguments_2_3_4'],
        'functions':event['functions'],
        'total_called_function_bytes':sum(x['bytes'] for x in event['functions']),
        'same_mask_group_address_IR_proof_reused':kind != 'Transport',
        'negative_ordinary_Default_witness':kind == 'Transport',
        'scope':'Original 4B guard exits before retained arithmetic. Exact prior optimized IR/source proof reused; new called machine bytes and allocator factory/engine provenance observed. No active multirow/numerical/benefit/frame claim.'}
out = {'structural_premise_observed':True, 'machine_policy_isolation_observed':True,
    'optimized_IR_unchanged':True, 'rows':rows,
    'actual_provider':{'path':str(driver.relative_to(R)),**pin(driver)},
    'intended_mixed_host':{'path':str(host.relative_to(R)),**pin(host)},
    'frozen_probe_is_intended_matched_host':False,
    'limits':'Factory selection and actual engine emission observed, not a full machine-pass trace or performance proof. No CPU/arithmetic/IR pipeline change inferred beyond exact saved-IR equality. All original C8/22/8 and changed-program runtime gates remain required.'}
target = W/'ordinary-code-review.json'; assert not target.exists()
target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'structural_premise_observed':True,'machine_policy_isolation_observed':True,
    'called_function_bytes':{k:v['total_called_function_bytes'] for k,v in rows.items()}}))
