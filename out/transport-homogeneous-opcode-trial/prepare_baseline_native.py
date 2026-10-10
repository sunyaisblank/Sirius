from pathlib import Path
import hashlib, importlib.util, json, subprocess

ROOT=Path.cwd(); WORK=Path(__file__).resolve().parent
revision=json.loads((WORK/'test-first-source.json').read_text())['source_revision']
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==revision
assert not subprocess.check_output(['git','status','--porcelain'])
output=WORK/'baseline-native-bindings.json'
assert not output.exists()
spec=importlib.util.spec_from_file_location('attestation',ROOT/'scripts/verify-attestation.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);g=v.load_build_gate_verifier()
stage=ROOT/'bin/windows-msvc'/('native-'+revision[:7])
bundle=ROOT/'attestations/native-build'/revision[:7]/'windows-build'
gate_path=stage/'generated/sirius/native_build_gate.json'
assert gate_path.read_bytes()==(bundle/'native_build_gate.json').read_bytes()
gate=g.validate_native_build_document(json.loads(gate_path.read_text()))
assert gate['source']=={'revision':revision,'clean':True}
v.verify_path(bundle/'windows-build.json',ROOT)
g.verify_recorded_files(gate,ROOT,stage)
assert len(v.copy_qualification_test_inputs(gate_path,stage))==20
pe=json.loads((WORK/'baseline-native-array-bindings.json').read_text())
assert pe['pass'] and pe['source_revision']==revision and len(pe['arrays'])==40 and not pe['changed_arrays']
files={}
def bind(path,expected=None):
    path=path.resolve();record={'bytes':path.stat().st_size,'sha256':sha(path)}
    if expected is not None: assert record=={key:expected[key] for key in ('bytes','sha256')},str(path)
    key=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
    assert key not in files or files[key]==record
    files[key]=record
for record in pe['products'].values():bind(ROOT/record['path'],record)
for path in stage.rglob('*'):
    assert not path.is_symlink()
    if path.is_file():bind(path)
for group in ('tested_artifacts','product_artifacts','test_input_artifacts'):
    for record in gate[group].values():bind((stage if record['root']=='build' else ROOT)/record['path'],record)
for path in subprocess.check_output(['git','ls-files','src/sirius/kernels','tests/support/retained_transport'],text=True).splitlines():bind(ROOT/path)
for path in ['src/sirius/backend/vulkan/vulkan_device.cpp','src/sirius/backend/retained_compute.cpp','src/sirius/backend/retained_compute.h','src/sirius/backend/retained_integrator.cpp','tests/backend/retained_compute_test.cpp','tests/backend/retained_dopri_test.cpp','scripts/verify-attestation.py','scripts/verify-build-gate.py','scripts/verify-operating-model.py','scripts/build-retained-kernels.py','.github/workflows/ci.yml','tests/labels/CTestLabels.cmake']:
    bind(ROOT/path)
for name in ['run_native_control.ps1','run_baseline_native.py','native_contracts.py','prepare_baseline_native.py','bind_baseline_native.py','baseline-terminal-audit.ps1','emergency_native_cleanup.ps1','test-first-source.json','baseline-build-and-arrays.json','baseline-retained_kernels.h','baseline-native-array-bindings.json','reconstruction-windows-'+revision[:7]+'.json','mixed-test-source-review.json','baseline-powershell-parser.json','baseline-native-apparatus-review.json','retention-gate-frozen.json','frozen-retention-comparator.py','test-first-publication.json','test-first-publication.txt']:
    bind(WORK/name)
bind(bundle/'windows-build.json')
provider={
'/mnt/c/Windows/System32/vulkan-1.dll':{'bytes':1719408,'sha256':'45ee2efc5c6d3986da8181775610f308122c752e55e703340d9733bbbc4b9f1a'},
'/mnt/c/Windows/System32/DriverStore/FileRepository/amdvlk.inf_amd64_914ba89eaaafdf60/amdvlk64.dll':{'bytes':103213064,'sha256':'5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650'},
'/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe':{'bytes':Path('/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe').stat().st_size,'sha256':'8bb6fa8c283b4d92120b1ef249a9b311b0f804d4cabbe9981159976c8be76a5e'}}
for path,record in provider.items():bind(Path(path),record)
result={'candidate_revision':revision,'baseline_revision':revision,'producers':{'baseline':{'source_revision':revision,'stage':str(stage.relative_to(ROOT)),'gate_sha256':sha(gate_path),'canonical_and_consumed_inputs_verified':20,'complete_recorded_files_verified':True}},'files':files,'file_count':len(files),'source_tree_clean':True,'qualification_claimed':False,'scope':'One new mixed-layer preservation baseline on exact accepted kernels; original official native products/all40 PE payloads/all20 consumed inputs/source/apparatus/provider pinned. Numerical execution and performance are pending.'}
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'files':len(files),'source_revision':revision,'dispatches':0}))
