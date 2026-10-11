"""Read one completed diagnostic against protected same-provider raw buffers."""
from pathlib import Path
import hashlib, json, subprocess
import software_contracts as contracts
WORK=Path(__file__).resolve().parent
ROOT=WORK.parents[1]
PRIOR=ROOT/'attestations/software-vulkan/dc714b9/guarded-upmultiply-trial/diagnostic'

def main():
    manifest_path=PRIOR.parent/'manifest.json'
    manifest=json.loads(manifest_path.read_text())
    for name in ('a1/gtest.xml','a1/owner.json','a1/inputs-before.json','baseline.json','expected-readback-layout.json','software_contracts.py'):
        matches=[r for r in manifest['payloads'] if r['path']=='diagnostic/'+name]
        assert len(matches)==1,name
        raw=(PRIOR/name).read_bytes();assert len(raw)==matches[0]['bytes'] and hashlib.sha256(raw).hexdigest()==matches[0]['sha256']
    original_owner=json.loads((PRIOR/'a1/owner.json').read_text());owner=json.loads((WORK/'candidate-timestamps/owner.json').read_text())
    for receipt in (owner,original_owner):
        assert receipt['passed'] and receipt['source_unchanged'] and receipt['whole_inputs_unchanged'] and receipt['child_exit']==0 and receipt['remaining_owned_processes']==[] and not receipt['cleanup_errors']
    current_inputs=json.loads((WORK/'candidate-timestamps/inputs-before.json').read_text());prior_inputs=json.loads((PRIOR/'a1/inputs-before.json').read_text())
    producer=json.loads((PRIOR/'baseline.json').read_text())
    original=(ROOT/'bin/linux-gcc/tests/backend/sirius_backend_tests').read_bytes()
    assert len(original)==producer['backend_bytes'] and hashlib.sha256(original).hexdigest()==producer['backend_sha256']
    assert subprocess.check_output(['git','rev-parse','HEAD^{tree}'],cwd=ROOT)==subprocess.check_output(['git','rev-parse',producer['source_revision']+'^{tree}'],cwd=ROOT)
    for name in ('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so','/usr/share/vulkan/icd.d/lvp_icd.json'):
        assert current_inputs[name]==prior_inputs[name],name
    contracts.ROOT=ROOT;contracts.WORK=WORK
    current=contracts.timestamp_result(contracts.xml_properties(WORK/'candidate-timestamps/gtest.xml','timestamps'),'candidate')
    contracts.WORK=PRIOR
    baseline=contracts.timestamp_result(contracts.xml_properties(PRIOR/'a1/gtest.xml','timestamps'),'baseline')
    assert current['readbacks']==baseline['readbacks'] and len(current['readbacks'])==32
    assert current['identity']==baseline['identity']
    assert {stage for stage in current['selected_modules'] if current['selected_modules'][stage]!=baseline['selected_modules'][stage]}=={'transport'}
    result={'pass':True,'scope':'Single current default diagnostic numerical case and protected historical same-provider full-buffer equality; no literal execution, performance gate or current-final qualification','current_owner_seconds':owner['wall_seconds'],'current_peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'current_input_seals':len(current_inputs),'full_readback_digests_exact':32,'ordered_available_markers':51,'selected_modules':current['selected_modules'],'provider_file_and_case_identity_exact':True,'historical_owner_revision':original_owner['source']['revision'],'historical_baseline_producer_revision':producer['source_revision'],'original_current_executable_byte_exact_historical_baseline_producer':True,'current_tracked_tree_exact_historical_baseline_producer':True,'historical_baseline_not_current_execution':True,'recognizer_result_observable_in_reserved_raw_word3':True,'canonical_marker0_preserves_original_bytes':True,'mismatch_marker1_deliberately_changes_raw_word3_not_executed_here':True,'physical_table_traffic_storage_and_cost_unresolved':True,'current_result':current}
    (WORK/'numerical-control-readback.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('pass','full_readback_digests_exact','ordered_available_markers','current_owner_seconds','current_peak_sampled_rss_kib','historical_owner_revision','historical_baseline_producer_revision')}))

if __name__=='__main__':
    main()
