"""Validate the existing CI export and reconstruct its exact ignored native stage."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, importlib.util, json, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HEAD = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert HEAD == '99a4a6d3a6bff69d9f20c995383a5f2e6b8cb901'
assert HEAD == subprocess.check_output(['git', 'rev-parse', '@{u}'], text=True).strip()
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
ARTIFACT = WORK / 'artifact'
STAGE = ROOT / 'bin/windows-msvc/native-99a4a6d-bound-select-trial'
assert not STAGE.exists()

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

preparation_review = json.loads((WORK/'native-preparation-review.json').read_text())
assert preparation_review['accepted'] and preparation_review['preparer_seal'] == seal(Path(__file__))
compiled_review = json.loads((WORK/'compiled-review.json').read_text())
assert compiled_review['accepted']
for name, key in (('product-inspection.json','product_inspection'), ('inspect_products.py','inspector'), ('candidate/retained_kernels.h','candidate_header')):
    assert seal(WORK/name) == compiled_review[key]
for name in ('source-review.json','product-inspection.json'):
    review = json.loads((WORK/name).read_text())
    assert review['accepted']
    for path, record in review['input_seals'].items():
        input_path = WORK/path if path == 'trial.patch' else ROOT/path
        assert seal(input_path) == record, path
controller_review = json.loads((WORK/'controller-review.json').read_text())
assert controller_review['accepted']
for name, record in controller_review['controller_seals'].items():
    assert seal(WORK/name) == record, name
analysis_review = json.loads((WORK/'analysis-review.json').read_text())
assert analysis_review['accepted']
for name, record in analysis_review['analysis_seals'].items():
    assert seal(WORK/name) == record, name
ci = json.loads((WORK/'ci-completed.json').read_text())
assert ci['databaseId'] == 38159862169 and ci['headSha'] == HEAD
assert ci['status'] == 'completed' and ci['conclusion'] == 'success'
assert {job['name'] for job in ci['jobs'] if job['conclusion'] == 'success'} == {'integration-no-render','integration-macos-no-render','integration-windows-no-render'}

def safe(root, relative):
    path = root / relative
    assert not Path(relative).is_absolute() and '..' not in Path(relative).parts
    assert path.resolve().is_relative_to(root.resolve())
    return path

document = json.loads((ARTIFACT / 'windows-build.json').read_text())
assert document['status'] == 'pass' and document['source_revision'] == HEAD
assert document['domains'] == ['windows-native-build']
assert document['platform'] == {'os': 'windows', 'wsl2': False}
assert document['claims']['compiled_topology_complete']
assert document['claims']['runtime_ready'] is False
payloads = {}
for name, record in document['artifacts'].items():
    path = safe(ARTIFACT, record['path'])
    assert seal(path) == {k: record[k] for k in ('bytes', 'sha256')}
    payloads.setdefault((record['bytes'], record['sha256']), path)

sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('native_gate_verifier', ROOT / 'scripts/verify-build-gate.py')
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)
gate = json.loads((ARTIFACT / 'native_build_gate.json').read_text())
verifier.validate_native_build_document(gate)
assert gate['source'] == {'revision': HEAD, 'clean': True}
STAGE.mkdir(parents=True)

def copy_record(record, destination):
    source = payloads[(record['bytes'], record['sha256'])]
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        assert destination.read_bytes() == source.read_bytes()
    else:
        shutil.copyfile(source, destination)
    assert seal(destination) == {k: record[k] for k in ('bytes', 'sha256')}

tested = {}
for name, record in gate['tested_artifacts'].items():
    assert record['root'] == 'build'
    tested[name] = safe(STAGE, record['path'])
    copy_record(record, tested[name])
for record in list(gate['product_artifacts'].values()) + list(gate['test_input_artifacts'].values()):
    if record['root'] == 'build':
        copy_record(record, safe(STAGE, record['path']))
    else:
        assert record['root'] == 'source'
        assert seal(safe(ROOT, record['path'])) == {k: record[k] for k in ('bytes', 'sha256')}

for name, record in gate['product_artifacts'].items():
    if name == 'sirius': continue
    relative = verifier.INSTALLED_PRODUCTS[name].removeprefix('share/sirius/')
    copy_record(record, safe(tested['sirius'].parent / 'resources', relative))
for consumer in ('sirius_backend_tests', 'sirius_render_tests'):
    resource = tested[consumer].parent / 'resources'
    for name, record in gate['test_input_artifacts'].items():
        copy_record(record, safe(resource, verifier.TEST_INPUT_PATHS[name]))
    for name in ('trace_spv', 'trace_fp32comp_spv', 'trace_fp64_spv'):
        copy_record(gate['product_artifacts'][name], safe(resource, verifier.INSTALLED_PRODUCTS[name].removeprefix('share/sirius/')))
verifier.verify_test_input_copies(tested, gate['product_artifacts'], gate['test_input_artifacts'])
gate_path = STAGE / 'generated/sirius/native_build_gate.json'
gate_path.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(ARTIFACT / 'native_build_gate.json', gate_path)

exe = tested['sirius_backend_tests'].read_bytes()
import re, struct
local = json.loads((WORK/'product-inspection.json').read_text())
assert local['accepted']
embedded = {}
for count, name, body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};', (WORK/'candidate/retained_kernels.h').read_text(),re.S):
    words = [int(word[:-1],0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)]
    assert len(words)==int(count)
    data = struct.pack('<'+'I'*len(words),*words)
    record = {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
    assert record==local['arrays'][name] and data in exe, name
    embedded[name] = {**record,'whole_payload_present_in_native_executable':True}
assert len(embedded)==40
baseline_stage = ROOT/'bin/windows-msvc/native-b81e061-executed-cost'
baseline_gate = json.loads((baseline_stage/'generated/sirius/native_build_gate.json').read_text())
assert baseline_gate['source']=={'clean':True,'revision':'b81e061fa081cdaed37eca3add950bd254bf296d'}
baseline = json.loads((WORK/'baseline-native-binding.json').read_text())
assert baseline['gate_source']==baseline_gate['source']
assert seal(baseline_stage/'generated/sirius/native_build_gate.json')==baseline['gate']
assert seal(baseline_stage/baseline['backend_executable']['path'])=={k:baseline['backend_executable'][k] for k in ('bytes','sha256')}
providers = json.loads((ROOT/'attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport/plan.json').read_text())['providers']
for provider in providers:
    assert seal(Path(provider['linux_path']))=={k:provider[k] for k in ('bytes','sha256')}
input_paths = [p for directory in (STAGE,baseline_stage) for p in directory.rglob('*') if p.is_file()]
input_paths += [WORK/name for name in ('prepare_native_stage.py','native-preparation-review.json','controller-review.json','native-cases.json','run_native.py','run_native_control.ps1','run_native_retirement.py','verify_native_retirement.ps1','emergency_native_cleanup.ps1','source-review.json','compiled-review.json','product-inspection.json','inspect_products.py','candidate/retained_kernels.h','baseline-native-binding.json','execution-preregistration.md','execution-preregistration-readback.json','ci-completed.json','artifact/windows-build.json','artifact/native_build_gate.json')]
input_paths += [ROOT/'attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport/plan.json', ROOT/'scripts/verify-build-gate.py']
input_paths += [WORK/name for name in ('compare_matched.py','analyze_controls.py','evaluate_retention.py','analysis-review.json')]
input_paths += [ROOT/name for name in ('src/sirius/backend/retained_compute.cpp','src/sirius/backend/vulkan/vulkan_device.cpp','tests/backend/retained_compute_test.cpp','tests/backend/retained_dopri_test.cpp','src/sirius/kernels/retained_pair.slang','scripts/build-retained-kernels.py')]
inputs = {str(p.relative_to(ROOT)):seal(p) for p in input_paths}
cases = json.loads((WORK/'native-cases.json').read_text())['cases']
actions = {name:{'filter':f,'producer':'candidate'} for name,f in cases.items()}
for label in ('a1','b1','b2','a2'):
    for case in ('factored','timestamps'):
        actions[label+'-'+case]={'filter':cases[case],'producer':'baseline' if label.startswith('a') else 'candidate'}
plan = {'source_revision':HEAD,'baseline_revision':'2053ef79c90752d37b0075b8137cfb9b571b6820','producers':{
 'candidate':{'source_revision':HEAD,'build_root':str(STAGE.relative_to(ROOT)),'executable_directory':str(tested['sirius_backend_tests'].parent.relative_to(ROOT)),'gate_sha256':seal(gate_path)['sha256']},
 'baseline':{'source_revision':baseline_gate['source']['revision'],'build_root':str(baseline_stage.relative_to(ROOT)),'executable_directory':str((baseline_stage/baseline['backend_executable']['path']).parent.relative_to(ROOT)),'gate_sha256':baseline['gate']['sha256']}},
 'input_seals':inputs,'providers':providers,'whole_candidate_arrays_in_native_executable':embedded,'cases':actions,
 'purpose':'Eight original candidate numerical/reference/orderedstream/refusal/recovery/rollback controls, then originalABBA sameprovider completebuffers and frozen44conditionwholecallbenefit gate. Original90s/4GiB nativechild guard,400s containedbridge; allscientific/device budgets/counts/governors unchanged. No fullqualification or adoptedoptimization.'}
assert HEAD==subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip() and not subprocess.check_output(['git','status','--porcelain'],text=True)
(WORK/'execution-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print(json.dumps({'source_revision':HEAD,'whole_native_arrays':len(embedded),'input_files':len(inputs),'native_backend':seal(tested['sirius_backend_tests']),'build_root':str(STAGE.relative_to(ROOT))}))
