from pathlib import Path
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET

work = Path(__file__).resolve().parent
root = Path.cwd()

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def git(*args):
    return subprocess.check_output(['git', *args], text=True).strip()

def birth(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None

source = json.loads((work/'source.json').read_text())
head = source['source_revision']
assert git('rev-parse', 'HEAD') == head and not git('status', '--porcelain')
build = json.loads((work/'build-and-arrays.json').read_text())
bindings = json.loads((work/'linux-readonly-payloads.json').read_text())
assert build['source_revision'] == bindings['source_revision'] == head
assert build['returncode'] == 0 and build['all_40_arrays_exact'] and bindings['pass_']
assert len(build['arrays']) == 40 and len(bindings['actual_consumers']) == 4
assert sha(root/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h') == build['whole_header_sha256'] == bindings['whole_header_sha256']
for relative, entry in bindings['actual_consumers'].items():
    data = (root/relative).read_bytes()
    assert len(data) == entry['bytes'] and hashlib.sha256(data).hexdigest() == entry['sha256']
    assert set(entry['arrays']) == set(build['arrays'])
    for name, record in entry['arrays'].items():
        assert {key:record[key] for key in ('bytes', 'sha256')} == build['arrays'][name]
        for occurrence in record['occurrences']:
            payload = data[occurrence['file_offset']:occurrence['file_offset']+record['bytes']]
            assert len(payload) == record['bytes'] and hashlib.sha256(payload).hexdigest() == record['sha256']

expected = {
 'preparation-model': ('RetainedComputeAdmission.SoftwareRendererPreparationPreservesPhysicalAccounting', 30),
 'renderer': ('VulkanRenderSession.ContinuationRendererPublishesOnlyCompleteFramesWithinActualBudget', 180),
}
complete = json.loads((work/'controls-complete.json').read_text())
assert complete['source_revision'] == head and complete['returncode'] == 0 and complete['existing_tests'] == 2
assert complete['new_test_registrations'] == 0 and not complete['full_qualification']
sequence = json.loads((work/'controls-sequence.json').read_text())
assert [item['name'] for item in sequence] == list(expected) and sequence == complete['results']
inputs = {}
births = {}
observations = []
for name, (case, limit) in expected.items():
    folder = work/'controls'/name
    owner = json.loads((folder/'owner.json').read_text())
    result = json.loads((folder/'result.json').read_text())
    assert result == next(item for item in sequence if item['name'] == name)
    assert owner['source_revision'] == head
    assert result['source_revision'] == head
    assert owner['case'] == case and owner['timeout_seconds'] == limit
    assert result['case'] == case
    assert owner['maximum_rss_mib'] == 4096 and owner['peak_sampled_rss_kib'] <= 4096*1024
    assert owner['elapsed_seconds'] < limit and owner['returncode'] == 0
    assert result['returncode'] == 0
    assert owner['stop_reason'] is None and owner['observer_error'] is None and not owner['cleanup_errors']
    assert not owner['gtest_environment_overrides']
    assert owner['owned_birth_absent'] and owner['owned_group_absent'] and owner['remaining_owned_group'] == []
    assert owner['source_and_inputs_unchanged']
    assert result['owner_sha256'] == sha(folder/'owner.json') and result['xml_sha256'] == sha(folder/'tests.xml')
    xml = ET.parse(folder/'tests.xml').getroot()
    assert [int(xml.get(key, '0')) for key in ['tests', 'failures', 'errors', 'disabled']] == [1, 0, 0, 0]
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('classname')+'.'+cases[0].get('name') == case
    assert cases[0].get('status') == 'run' and cases[0].get('result') == 'completed' and not cases[0].findall('skipped')
    assert not cases[0].findall('failure')
    if name not in ['preparation-model', 'wire']:
        assert any('libvulkan_lvp.so' in path for path in owner['sampled_loaded_provider_paths'])
        assert any('libvulkan.so' in path for path in owner['sampled_loaded_provider_paths'])
    for path, identity in owner['whole_input_seals'].items():
        if path in inputs: assert inputs[path] == identity
        inputs[path] = identity
    births[(owner['controller_pid'], owner['controller_start_ticks'])] = 'controller'
    for entry in owner['saved_birth_checks']:
        assert entry['absent']
        births[(entry['pid'], entry['saved_start_ticks'])] = 'test/group'
    observations.append({'name':name, 'case':case, 'elapsed_seconds':owner['elapsed_seconds'],
                         'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'], 'zero_failure_error_skip':True})
for path, identity in inputs.items():
    artifact = Path(path)
    assert artifact.stat().st_size == identity['bytes'] and sha(artifact) == identity['sha256']
for command in build['commands']:
    assert command['returncode'] == 0
    births[(command['pid'], command['start_ticks'])] = 'configure/build'
birth_checks = [{'pid':pid, 'start_ticks':saved, 'role':role, 'absent':birth(pid) != saved}
                for (pid, saved), role in sorted(births.items())]
assert all(item['absent'] for item in birth_checks)
assert git('rev-parse', 'HEAD') == head and not git('status', '--porcelain')
receipt = {'source_revision':head, 'source_tree':source['source_tree'], 'verification_pass':True, 'batch_pass':True, 'passed_cases':2, 'failed_cases':0,
           'existing_cases':2, 'new_registrations':0, 'observations':observations,
           'selected_whole_inputs_postverified':len(inputs), 'birth_checks':birth_checks,
           'actual_elf_consumers':4, 'complete_array_bindings':160, 'all_40_arrays_exact':True,
           'verifier_sha256':sha(Path(__file__)), 'full_qualification':False,
           'scope':'Postprocessing of exactly two passing preparation-count correction controls and actual whole input/array bytes. No new numerical execution, speed benefit, native hardware, required failed full frame or complete scientific/release qualification claim.'}
(work/'verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['source_revision', 'verification_pass', 'batch_pass', 'passed_cases', 'failed_cases', 'existing_cases', 'selected_whole_inputs_postverified', 'complete_array_bindings']}))
