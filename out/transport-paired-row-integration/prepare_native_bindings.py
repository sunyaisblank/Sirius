"""Bind each exact native producer to its own original authority and controllers."""
from pathlib import Path
import hashlib, importlib.util, json, subprocess, sys

ROOT=Path.cwd(); WORK=Path(__file__).resolve().parent
mode,=sys.argv[1:]; assert mode in ('baseline','execution')
source=json.loads((WORK/'source.json').read_text()); head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def document(path):return json.loads(path.read_text(encoding='utf-8-sig'))
files={}
def bind(path,expected=None):
    path=path.resolve();assert path.is_file()
    record={'bytes':path.stat().st_size,'sha256':sha(path)}
    if expected is not None:assert record=={k:expected[k] for k in ('bytes','sha256')},str(path)
    key=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
    assert key not in files or files[key]==record;files[key]=record
spec=importlib.util.spec_from_file_location('attestation',ROOT/'scripts/verify-attestation.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);g=v.load_build_gate_verifier()
producers={}
for label,revision in [('baseline',source['baseline_revision']),('candidate',head)]:
    stage=ROOT/'bin/windows-msvc'/('native-'+revision[:7])
    bundle=ROOT/'attestations/native-build'/revision[:7]/'windows-build'
    receipt=document(WORK/('reconstruction-windows-'+revision[:7]+'.json'))
    assert receipt['source_revision']==revision and receipt['live_source_revision']==head
    view=ROOT/receipt['recorded_source_view'];assert view==stage/'recorded-source'
    gate_path=stage/'generated/sirius/native_build_gate.json'
    assert gate_path.read_bytes()==(bundle/'native_build_gate.json').read_bytes()
    assert sha(gate_path)==receipt['gate_sha256']
    gate=g.validate_native_build_document(document(gate_path))
    assert gate['source']=={'revision':revision,'clean':True}
    model=subprocess.check_output(['git','show',revision+':tests/operating_model.json'])
    v.verify_document_against_authority(document(bundle/'windows-build.json'),bundle/'windows-build.json',revision,hashlib.sha256(model).hexdigest())
    for collection in ('tested_artifacts','product_artifacts','test_input_artifacts'):
        for record in gate[collection].values():
            path=(stage if record['root']=='build' else view)/record['path']
            if record['root']=='source':
                assert path.read_bytes()==subprocess.check_output(['git','show',revision+':'+record['path']])
            bind(path,record)
    g.verify_recorded_files(gate,view,stage)
    assert len(v.copy_qualification_test_inputs(gate_path,stage))==20
    pe=document(WORK/(label+'-native-array-bindings.json'))
    assert pe['pass'] and pe['source_revision']==revision and len(pe['arrays'])==40
    assert pe['changed_arrays']==([] if label=='baseline' else ['kTransportShader'])
    for record in pe['products'].values():bind(ROOT/record['path'],record)
    exported=document(WORK/(label+'-windows-export-download.json'))
    assert exported['source_revision']==revision and len(exported['extracted_files'])==40
    assert exported['artifact']['digest']=='sha256:'+exported['sha256']
    assert {q['path'] for q in exported['extracted_files']}=={str(p.relative_to(bundle)) for p in bundle.rglob('*') if p.is_file()}
    for record in exported['extracted_files']:bind(bundle/record['path'],record)
    for path in stage.rglob('*'):
        assert not path.is_symlink()
        if path.is_file():bind(path)
    producers[label]={'source_revision':revision,'stage':str(stage.relative_to(ROOT)),
                      'gate_sha256':sha(gate_path),'recorded_source_view':str(view.relative_to(ROOT)),
                      'canonical_and_consumed_inputs_verified':20,'complete_recorded_files_verified':True}
    for name in [label+'-native-array-bindings.json',label+'-windows-export-download.json',
                 label+'-native-ci-artifacts.json',label+'-native-export-ci-readback.json',
                 'reconstruction-windows-'+revision[:7]+'.json']:bind(WORK/name)
for path in subprocess.check_output(['git','ls-files','src/sirius/kernels','tests/support/retained_transport'],text=True).splitlines():bind(ROOT/path)
for path in [*source['changed_files'],'src/sirius/backend/vulkan/vulkan_device.cpp',
             'src/sirius/backend/retained_integrator.cpp','src/sirius/backend/retained_integrator.h',
             'src/sirius/backend/retained_trace_executor.cpp','src/sirius/backend/retained_trace_executor.h',
             'tests/backend/retained_dopri_test.cpp','scripts/verify-attestation.py',
             'scripts/verify-build-gate.py','.github/workflows/ci.yml','CMakePresets.json']:bind(ROOT/path)
for name in ['source.json','build-and-arrays.json','linux-readonly-payloads.json','controls-complete.json',
             'independent-software-controls-review.json','native-export-apparatus-review.json',
             'candidate-native-apparatus-review.json','native-powershell-parser-reference.json',
             'original-timestamp-reference.json','retention-gate-frozen.json','independent-freeze-review.json',
             'run_native.py','run_native_sequence.py','run_terminal_audit.py','native_contracts.py',
             'run_native_control.ps1','baseline-terminal-audit.ps1','emergency_native_cleanup.ps1',
             'frozen-retention-comparator.py','download_exports.py','reconstruct_native.py','bind_native.py',
             'prepare_native_bindings.py','candidate-native-preregistration.txt','candidate-native-preregistration.json',
             'controls/timestamps/tests.xml']:bind(WORK/name)
assert document(WORK/'controls-complete.json')['returncode']==0
assert document(WORK/'independent-software-controls-review.json')['pass_']
assert document(WORK/'native-export-apparatus-review.json')['pass_']
assert document(WORK/'candidate-native-apparatus-review.json')['pass']
assert document(WORK/'independent-freeze-review.json')['pass_']
assert document(WORK/'candidate-native-preregistration.json')['readback']['body']==(WORK/'candidate-native-preregistration.txt').read_text()
freeze=document(WORK/'retention-gate-frozen.json')
assert sha(WORK/'frozen-retention-comparator.py')==freeze['comparator_sha256']
assert freeze['candidate_revision']==head and freeze['baseline_revision']==source['baseline_revision']
provider={
 '/mnt/c/Windows/System32/vulkan-1.dll':{'bytes':1719408,'sha256':'45ee2efc5c6d3986da8181775610f308122c752e55e703340d9733bbbc4b9f1a'},
 '/mnt/c/Windows/System32/DriverStore/FileRepository/amdvlk.inf_amd64_914ba89eaaafdf60/amdvlk64.dll':{'bytes':103213064,'sha256':'5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650'},
 '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe':{'bytes':495616,'sha256':'8bb6fa8c283b4d92120b1ef249a9b311b0f804d4cabbe9981159976c8be76a5e'}}
for path,record in provider.items():bind(Path(path),record)
if mode=='execution':
    accepted=WORK/'terminal-baseline-baseline-mixed-accepted.json'
    accepted_record=document(accepted);assert accepted_record['pass'];bind(accepted)
    for artifact in accepted_record['record']['artifacts']:bind(ROOT/artifact['path'],artifact)
    terminal=ROOT/accepted_record['terminal_receipt'];assert sha(terminal)==accepted_record['terminal_receipt_sha256'];bind(terminal)
    terminal_bridge=WORK/'terminal-baseline-baseline-mixed-bridge.json';assert sha(terminal_bridge)==accepted_record['terminal_bridge_sha256'];bind(terminal_bridge)
    baseline=document(WORK/'baseline-native-bindings.json')
    assert baseline['producers']==producers
    for relative,record in baseline['files'].items():bind(Path(relative) if relative.startswith('/') else ROOT/relative,record)
    bind(WORK/'baseline-native-bindings.json')
output=WORK/('baseline-native-bindings.json' if mode=='baseline' else 'execution-bindings.json')
assert not output.exists()
record={'candidate_revision':head,'baseline_revision':source['baseline_revision'],'producers':producers,
        'files':files,'file_count':len(files),'source_tree_clean':True,'qualification_claimed':False,
        'scope':'Exact own-commit official native build gates/products/resources/20consumedinputs; actual40-array wholePEbindings, relevant source/controllers/frozencriteria/preregistration and selected providers. No native outcome/adoption/frame/fullqualification claim.'}
output.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'mode':mode,'files':len(files),'producers':producers,'dispatches':0}))
