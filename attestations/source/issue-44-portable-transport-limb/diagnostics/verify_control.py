"""Read the unchanged integrated contracts and exact complete-buffer metadata."""
from pathlib import Path
import hashlib
import importlib.util
import json
import resource
import sys

sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024*1024,64*1024*1024))
root=Path.cwd()
work=Path(__file__).resolve().parent
metadata=work/'control-contracts'
spec=importlib.util.spec_from_file_location('unchanged_original_contract_reader',metadata/'software_contracts.py')
reader=importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)
assert reader.ROOT==root and reader.WORK==metadata
folder=work/'integrated-control-8'
owner=json.loads((folder/'owner.json').read_text())
assert owner['passed'] and owner['status']=='terminal' and owner['whole_inputs_unchanged']
assert owner['source_unchanged'] and owner['remaining_owned_processes']==[] and not owner['cleanup_errors']
props=reader.xml_properties(folder/'gtest.xml','timestamps')
current=reader.timestamp_result(props,'candidate')
old=root/'attestations/software-vulkan/41fe5f2/portable-denominator-preparation/files/out/portable-denominator-preparation/benefit/a1/gtest.xml'
baseline=reader.timestamp_result(reader.xml_properties(old,'timestamps'),'baseline')
assert current['identity']==baseline['identity']
assert current['readbacks']==baseline['readbacks'],'Current candidate changed an original complete-buffer digest'
assert len(current['readbacks'])==32 and len(current['executed_modules'])==54
assert len(current['selected_modules'])==7
report={'passed':True,'gtest_tests':1,'failures':0,'errors':0,'skipped':0,
    'original32complete_buffer_digests_unchanged':True,'ordered_markers':51,
    'actual_command_module_joins':54,'selected_module_joins':7,
    'reader_sha256':hashlib.sha256((metadata/'software_contracts.py').read_bytes()).hexdigest(),
    'current_xml_sha256':hashlib.sha256((folder/'gtest.xml').read_bytes()).hexdigest(),
    'historical_baseline_xml_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
    'identity':current['identity'],'readbacks':current['readbacks'],
    'actual_executed_modules':current['executed_modules'],
    'scope':'Original integratedcontrol8 numerical/source/stream/refusal/recovery plus actualcurrentmodulejoins. '
        'Historical baseline XML is digest/layout identity evidence only, never a current performance target. '
        'Sole frozen ABBA/remaining7controls/fullscientific/frame/releasequalification remain pending.'}
output=work/'integrated-control-check.json'
assert not output.exists()
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('identity','readbacks','actual_executed_modules')}))
