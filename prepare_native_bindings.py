"""Bind exact own-commit native producers and the reviewed default apparatus."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import hashlib, importlib.util, json, subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
label,=sys.argv[1:];assert label in ('baseline','candidate')
def document(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=document(WORK/'source.json');head=source['source_revision'];assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head;assert not subprocess.check_output(['git','status','--porcelain'])
assert document(WORK/(label+'-native-apparatus-review.json'))['pass']
assert document(WORK/'independent-build-admission-review.json')['pass']
assert document(WORK/'candidate-native-preregistration.json')['readback']==document(WORK/'issue-readback.json')
assert document(WORK/'issue-readback.json')['body']==(WORK/'candidate-native-preregistration.txt').read_text()
assert document(WORK/'retention-gate-frozen.json')['checks']==46
assert sha(WORK/'compare_native.py')==document(WORK/'retention-gate-frozen.json')['comparator_sha256']=='7e13981d4fb3e566c5e5ae404946a899530096d0d3c6f7b09c664fc5f8a5f686'
spec=importlib.util.spec_from_file_location('attestation',ROOT/'scripts/verify-attestation.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);g=v.load_build_gate_verifier()
files={};producers={}
def bind(path,expected=None):
 path=path.resolve();assert path.is_file() and not path.is_symlink();value={'bytes':path.stat().st_size,'sha256':sha(path)}
 if expected is not None:assert value=={k:expected[k] for k in ('bytes','sha256')},str(path)
 key=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path);assert key not in files or files[key]==value;files[key]=value
requested=[('baseline',source['baseline_revision'])]+([] if label=='baseline' else [('candidate',head)])
for producer,revision in requested:
 stage=ROOT/'bin/windows-msvc'/('native-'+revision[:7]);bundle=ROOT/'attestations/native-build'/revision[:7]/'windows-build';receipt=document(WORK/('reconstruction-windows-'+revision[:7]+'.json'))
 assert receipt['source_revision']==revision and receipt['live_source_revision']==head and receipt['canonical_and_consumed_inputs_checked']==20
 view=ROOT/receipt['recorded_source_view'];assert view==stage/'recorded-source';gate_path=stage/'generated/sirius/native_build_gate.json';assert gate_path.read_bytes()==(bundle/'native_build_gate.json').read_bytes() and sha(gate_path)==receipt['gate_sha256']
 gate=g.validate_native_build_document(document(gate_path));assert gate['source']=={'revision':revision,'clean':True}
 model=subprocess.check_output(['git','show',revision+':tests/operating_model.json']);v.verify_document_against_authority(document(bundle/'windows-build.json'),bundle/'windows-build.json',revision,hashlib.sha256(model).hexdigest())
 for collection in ('tested_artifacts','product_artifacts','test_input_artifacts'):
  for record in gate[collection].values():
   path=(stage if record['root']=='build' else view)/record['path']
   if record['root']=='source':assert path.read_bytes()==subprocess.check_output(['git','show',revision+':'+record['path']])
   bind(path,record)
 g.verify_recorded_files(gate,view,stage);assert len(v.copy_qualification_test_inputs(gate_path,stage))==20
 pe=document(WORK/(producer+'-native-array-bindings.json'));assert pe['pass'] and pe['source_revision']==revision and len(pe['arrays'])==(40 if producer=='baseline' else 42) and pe['original40_unchanged']
 assert pe['whole_payload_joins']==(120 if producer=='baseline' else 126)
 for record in pe['products'].values():bind(ROOT/record['path'],record)
 for directory in (stage,bundle):
  for path in sorted(directory.rglob('*')):
   assert not path.is_symlink()
   if path.is_file():bind(path)
 for name in (producer+'-native-array-bindings.json','reconstruction-windows-'+revision[:7]+'.json'):bind(WORK/name)
 producers[producer]={'source_revision':revision,'stage':str(stage.relative_to(ROOT)),'gate_sha256':sha(gate_path),'recorded_source_view':str(view.relative_to(ROOT)),'canonical_and_consumed_inputs_verified':20,'complete_recorded_files_verified':True}
if label=='candidate':
 from native_contracts import read_case
 bind(WORK/'baseline-native-bindings.json')
 baseline_bindings=document(WORK/'baseline-native-bindings.json');observed=read_case('baseline','baseline','timestamps',baseline_bindings);accepted=document(WORK/'terminal-baseline-baseline-timestamps-accepted.json')
 assert accepted['pass'] and accepted['record']['artifacts']==observed['artifacts']
 admission=document(WORK/'baseline-native-admission-accepted.json');assert admission['pass'] and admission['observed']==observed and admission['terminal_acceptance']==accepted
 for path in sorted((WORK/'native/baseline/baseline/timestamps').iterdir()):
  if path.is_file():bind(path)
for path in subprocess.check_output(['git','ls-files'],text=True).splitlines():bind(ROOT/path)
# All currently completed root receipts/controllers are immutable; later native outputs are separate paths.
for path in sorted(WORK.iterdir()):
 if path.is_file() and path.name not in ('execution-bindings.json','baseline-native-bindings.json'):bind(path)
for path,record in document(WORK/'provider-current.json')['files'].items():bind(Path(path),record)
layout=document(WORK/'expected-readback-layout.json');bind(ROOT/layout['source'],{'bytes':layout['source_bytes'],'sha256':layout['source_sha256']})
freeze=document(WORK/'retention-gate-frozen.json');bind(ROOT/freeze['original_gate_manifest'],{'bytes':(ROOT/freeze['original_gate_manifest']).stat().st_size,'sha256':freeze['original_gate_manifest_sha256']});bind(ROOT/freeze['comparator_reference'],{'bytes':(ROOT/freeze['comparator_reference']).stat().st_size,'sha256':freeze['comparator_sha256']})
output=WORK/('baseline-native-bindings.json' if label=='baseline' else 'execution-bindings.json');assert not output.exists()
output.write_text(json.dumps({'candidate_revision':head,'baseline_revision':source['baseline_revision'],'producers':producers,'files':files,'file_count':len(files),'source_tree_clean':True,'qualification_claimed':False,'fresh_baseline_logical_fma_required_before_candidate':True,'scope':'Own-commit official source/product/input authority and complete baseline40/candidate42 PE payloads; original40 unchanged; exact default46 comparator and reviewed controllers/provider bytes. Native result/benefit/full qualification unclaimed.'},indent=2)+'\n');print(json.dumps({'label':label,'files':len(files),'producers':producers,'dispatches':0}))
