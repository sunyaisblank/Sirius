import pathlib,json,hashlib,importlib.util,subprocess,shutil,sys
root=pathlib.Path.cwd();platform,label=sys.argv[1:];work=pathlib.Path(__file__).resolve().parent
source=json.loads((work/'source.json').read_text());assert label in ['baseline','candidate'];revision=source['baseline_revision'] if label=='baseline' else source['source_revision'];assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==source['source_revision'];assert not subprocess.check_output(['git','status','--porcelain'])
assert platform in ['windows','macos'] and len(revision)==40
bundle=root/'attestations/native-build'/revision[:7]/(platform+'-build')
stage=root/'bin'/('windows-msvc' if platform=='windows' else 'macos-clang')/('native-'+revision[:7])
spec=importlib.util.spec_from_file_location('sirius_attestation',root/'scripts/verify-attestation.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
g=v.load_build_gate_verifier();data=json.loads((bundle/(platform+'-build.json')).read_text());gate=g.validate_native_build_document(json.loads((bundle/'native_build_gate.json').read_text()))
assert gate['source']=={'revision':revision,'clean':True} and data['source_revision']==revision
assert subprocess.check_output(['git','cat-file','-t',revision],text=True).strip()=='commit'
model=subprocess.check_output(['git','show',revision+':tests/operating_model.json']);assert label=='baseline' or model==(root/'tests/operating_model.json').read_bytes()
head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
if head==revision:
 v.verify_path(bundle/(platform+'-build.json'),root)
else:
 # A preserved baseline keeps its own committed model and source authority.
 v.verify_document_against_authority(data,bundle/(platform+'-build.json'),revision,hashlib.sha256(model).hexdigest())
assert not stage.exists();stage.mkdir(parents=True)
source_view=stage/'recorded-source';source_differences=[];source_records=[]
def put(src,relative):
 rel=pathlib.Path(relative);assert not rel.is_absolute() and '..' not in rel.parts
 dst=stage/rel;dst.parent.mkdir(parents=True,exist_ok=True)
 if dst.exists():assert dst.read_bytes()==src.read_bytes()
 else:shutil.copy2(src,dst)
 return dst
def bound(src,record):
 b=src.read_bytes();assert len(b)==record['bytes'] and hashlib.sha256(b).hexdigest()==record['sha256']
 if record['root']=='build':return put(src,record['path'])
 assert record['root']=='source'
 canonical=subprocess.check_output(['git','show',revision+':'+record['path']])
 assert canonical==b
 if canonical!=(root/record['path']).read_bytes():source_differences.append(record['path'])
 source_records.append(record)
 return put(src,pathlib.Path('recorded-source')/record['path'])
executables={}
for name,record in gate['tested_artifacts'].items():
 src=bundle/('qualification-sirius.bin' if name=='sirius' else v.QUALIFICATION_TEST_EVIDENCE[name]);executables[name]=bound(src,record)
for name,record in gate['product_artifacts'].items():
 src=bundle/('qualification-sirius.bin' if name=='sirius' else ('alignment_receipt.json' if name=='alignment_receipt' else v.QUALIFICATION_PRODUCT_EVIDENCE[name]));bound(src,record)
for name,record in gate['test_input_artifacts'].items():bound(bundle/v.QUALIFICATION_TEST_INPUT_EVIDENCE[name],record)
for key,name in [('junit','native_build_gate_junit.xml'),('log','native-build-gate-log')]:bound(bundle/name,gate['ctest'][key])
put(bundle/'native_build_gate.json','generated/sirius/native_build_gate.json')
for consumer in ['sirius','sirius_backend_tests','sirius_render_tests']:
 parent=executables[consumer].parent.relative_to(stage)
 for name,rel in v.QUALIFICATION_RUNTIME_RESOURCE_PATHS.items():put(bundle/v.QUALIFICATION_PRODUCT_EVIDENCE[name],parent/'resources'/rel)
 put(bundle/'alignment_receipt.json',parent/'resources/model/alignment_receipt.json')
 if consumer!='sirius':
  for name,rel in g.TEST_INPUT_PATHS.items():put(bundle/v.QUALIFICATION_TEST_INPUT_EVIDENCE[name],parent/'resources'/rel)
assert source_differences==[], source_differences
g.verify_recorded_files(gate,source_view,stage)
checked=v.copy_qualification_test_inputs(stage/'generated/sirius/native_build_gate.json',stage);assert len(checked)==20
receipt={'source_revision':revision,'live_source_revision':head,'original_ci_export':str(bundle.relative_to(root)),'whole_recorded_files_verified':True,'source_root_records_match_committed_bytes':True,'recorded_source_view':str(source_view.relative_to(root)),'source_records':source_records,'intentional_baseline_local_differences':source_differences,'canonical_and_consumed_inputs_checked':20,'source_eol_exceptions':[],'full_native_runtime_qualification_claimed':False,'ctest_registration_reconstructed':False,'executables':{n:{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for n,p in executables.items()},'gate_sha256':hashlib.sha256((bundle/'native_build_gate.json').read_bytes()).hexdigest()}
(work/('reconstruction-'+platform+'-'+revision[:7]+'.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'platform':platform,'source_revision':revision,'whole_records_verified':True,'source_eol_exceptions':0,'consumed_inputs':20}))
