"""Bind exact native producers and the existing finite native apparatus."""
from pathlib import Path
import hashlib, importlib.util, json, subprocess
ROOT=Path.cwd(); WORK=Path(__file__).resolve().parent
def document(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
source=document(WORK/'source.json'); head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
assert document(WORK/'candidate-native-apparatus-review.json')['pass']
spec=importlib.util.spec_from_file_location('attestation',ROOT/'scripts/verify-attestation.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);g=v.load_build_gate_verifier()
files={};producers={}
def bind(path,expected=None):
    path=path.resolve(); assert path.is_file() and not path.is_symlink()
    value={'bytes':path.stat().st_size,'sha256':sha(path)}
    if expected is not None: assert value=={k:expected[k] for k in ('bytes','sha256')},str(path)
    key=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
    assert key not in files or files[key]==value;files[key]=value
for label,revision in [('baseline',source['baseline_revision']),('candidate',head)]:
    stage=ROOT/'bin/windows-msvc'/('native-'+revision[:7])
    bundle=ROOT/'attestations/native-build'/revision[:7]/'windows-build'
    receipt=document(WORK/('reconstruction-windows-'+revision[:7]+'.json'))
    assert receipt['source_revision']==revision and receipt['live_source_revision']==head
    view=ROOT/receipt['recorded_source_view'];assert view==stage/'recorded-source'
    gate_path=stage/'generated/sirius/native_build_gate.json'
    assert gate_path.read_bytes()==(bundle/'native_build_gate.json').read_bytes()
    assert sha(gate_path)==receipt['gate_sha256']
    gate=g.validate_native_build_document(document(gate_path));assert gate['source']=={'revision':revision,'clean':True}
    model=subprocess.check_output(['git','show',revision+':tests/operating_model.json'])
    v.verify_document_against_authority(document(bundle/'windows-build.json'),bundle/'windows-build.json',revision,hashlib.sha256(model).hexdigest())
    for collection in ('tested_artifacts','product_artifacts','test_input_artifacts'):
        for record in gate[collection].values():
            path=(stage if record['root']=='build' else view)/record['path']
            if record['root']=='source': assert path.read_bytes()==subprocess.check_output(['git','show',revision+':'+record['path']])
            bind(path,record)
    g.verify_recorded_files(gate,view,stage)
    assert len(v.copy_qualification_test_inputs(gate_path,stage))==20
    pe=document(WORK/(label+'-native-array-bindings.json'))
    assert pe['pass'] and pe['source_revision']==revision and len(pe['arrays'])==(40 if label=='baseline' else 42) and pe['original40_unchanged']
    for record in pe['products'].values():bind(ROOT/record['path'],record)
    for directory in (stage,bundle):
        for path in directory.rglob('*'):
            assert not path.is_symlink()
            if path.is_file():bind(path)
    for name in (label+'-native-array-bindings.json','reconstruction-windows-'+revision[:7]+'.json'):bind(WORK/name)
    producers[label]={'source_revision':revision,'stage':str(stage.relative_to(ROOT)),'gate_sha256':sha(gate_path),'recorded_source_view':str(view.relative_to(ROOT)),'canonical_and_consumed_inputs_verified':20,'complete_recorded_files_verified':True}
for path in subprocess.check_output(['git','ls-files'],text=True).splitlines():bind(ROOT/path)
for path in sorted(WORK.iterdir()):
    if path.is_file():bind(path)
bind(WORK/'software-wide-observer/gtest.xml')
bind(WORK/'software-wide-observer/accepted.json')
for path,expected in {
 '/mnt/c/Windows/System32/vulkan-1.dll':{'bytes':1719408,'sha256':'45ee2efc5c6d3986da8181775610f308122c752e55e703340d9733bbbc4b9f1a'},
 '/mnt/c/Windows/System32/DriverStore/FileRepository/amdvlk.inf_amd64_914ba89eaaafdf60/amdvlk64.dll':{'bytes':103213064,'sha256':'5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650'},
 '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe':{'bytes':495616,'sha256':'8bb6fa8c283b4d92120b1ef249a9b311b0f804d4cabbe9981159976c8be76a5e'}
}.items():bind(Path(path),expected)
freeze=document(WORK/'retention-gate-frozen.json')
assert freeze['source_revision']==head and freeze['baseline_revision']==source['baseline_revision']
assert freeze['comparator_sha256']==sha(WORK/'compare_native.py')
issue=document(WORK/'issue-readback.json')
assert issue['body']==(WORK/'candidate-native-preregistration.txt').read_text()
assert document(WORK/'candidate-native-preregistration.json')['readback']==issue
output=WORK/'execution-bindings.json';assert not output.exists()
output.write_text(json.dumps({'candidate_revision':head,'baseline_revision':source['baseline_revision'],'producers':producers,'files':files,'file_count':len(files),'source_tree_clean':True,'qualification_claimed':False,'scope':'Own-commit official product/source/input authority; complete baseline40/candidate42 PEpayloads, alloriginal40 unchanged; frozen source/controllers and actual selectedprovider modules. No native result, latency or fullqualification claim.'},indent=2)+'\n')
print(json.dumps({'files':len(files),'producers':producers,'dispatches':0}))
