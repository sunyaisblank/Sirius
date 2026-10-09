import pathlib, json, shutil, importlib.util, hashlib, subprocess
root=pathlib.Path.cwd(); revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert revision=='6e7ba9347a6eb15ccc607b85917edb7d9be5c34b'
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
bundle=root/'attestations/native-build/6e7ba93/windows-build'; stage=root/'bin/windows-msvc/native-6e7ba93'
subprocess.run(['python3',str(root/'scripts/verify-attestation.py'),str(bundle/'windows-build.json'),'--source-root',str(root)],check=True)
spec=importlib.util.spec_from_file_location('sirius_attestation',root/'scripts/verify-attestation.py'); v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
g=v.load_build_gate_verifier(); gate=json.loads((bundle/'native_build_gate.json').read_text()); assert gate['source']=={'revision':revision,'clean':True}
stage.mkdir(parents=True,exist_ok=True)
source_eol=[]
def put(src, relative):
 p=pathlib.Path(relative);assert not p.is_absolute() and '..' not in p.parts
 dst=stage/p;dst.parent.mkdir(parents=True,exist_ok=True)
 if dst.exists(): assert src.read_bytes()==dst.read_bytes()
 else:shutil.copy2(src,dst)
 return dst
def bound_copy(src,record):
 b=src.read_bytes();assert len(b)==record['bytes'] and hashlib.sha256(b).hexdigest()==record['sha256']
 if record['root']=='build':return put(src,record['path'])
 assert record['root']=='source'
 local=(root/record['path']).read_bytes()
 if local!=b:
  assert record['path'] in ['src/sirius/app/viewer/shaders/RDSD003A.frag','src/sirius/app/viewer/shaders/RDSD003A.vert']
  committed=subprocess.check_output(['git','show',revision+':'+record['path']])
  assert local==committed and b.replace(b'\r\n',b'\n')==committed and b'\r\n' in b
  source_eol.append({'path':record['path'],'original_export_sha256':record['sha256'],'local_git_blob_sha256':hashlib.sha256(committed).hexdigest(),'scope':'Exact exported CRLF bytes retained; canonical local LF source agrees with committed Git content. Source-root byte-exact gate reconstruction is not claimed.'})
 return root/record['path']
executables={}
for name,record in gate['tested_artifacts'].items():
 src=bundle/('qualification-sirius.bin' if name=='sirius' else v.QUALIFICATION_TEST_EVIDENCE[name]);executables[name]=bound_copy(src,record)
for name,record in gate['product_artifacts'].items():
 src=bundle/('qualification-sirius.bin' if name=='sirius' else ('alignment_receipt.json' if name=='alignment_receipt' else v.QUALIFICATION_PRODUCT_EVIDENCE[name]));bound_copy(src,record)
for name,record in gate['test_input_artifacts'].items():bound_copy(bundle/v.QUALIFICATION_TEST_INPUT_EVIDENCE[name],record)
for key in ['junit','log']:
 src=bundle/(gate['ctest'][key]['path'].split('/')[-1] if key=='junit' else 'native-build-gate-log');bound_copy(src,gate['ctest'][key])
put(bundle/'native_build_gate.json','generated/sirius/native_build_gate.json')
for name,relative in v.QUALIFICATION_RUNTIME_RESOURCE_PATHS.items():
 src=bundle/v.QUALIFICATION_PRODUCT_EVIDENCE[name];put(src,executables['sirius'].parent.relative_to(stage)/'resources'/relative)
put(bundle/'alignment_receipt.json',executables['sirius'].parent.relative_to(stage)/'resources/model/alignment_receipt.json')
for consumer in ['sirius_backend_tests','sirius_render_tests']:
 for name,relative in g.TEST_INPUT_PATHS.items():
  src=bundle/v.QUALIFICATION_TEST_INPUT_EVIDENCE[name];put(src,executables[consumer].parent.relative_to(stage)/'resources'/relative)
 # Runtime products travel with each consumer too, as the existing stage-runtime helper does.
 for name,relative in v.QUALIFICATION_RUNTIME_RESOURCE_PATHS.items():put(bundle/v.QUALIFICATION_PRODUCT_EVIDENCE[name],executables[consumer].parent.relative_to(stage)/'resources'/relative)
 put(bundle/'alignment_receipt.json',executables[consumer].parent.relative_to(stage)/'resources/model/alignment_receipt.json')
# The original complete bundle verifies against the source; replay every staged
# build-root record. Two CI source-root viewer records legitimately carry CRLF,
# so whole verify_recorded_files against the Linux checkout is not claimed.
for collection in ['tested_artifacts','product_artifacts','test_input_artifacts']:
 for record in gate[collection].values():
  if record['root']=='build':g.resolve_record(record,root,stage)
for key in ['junit','log']:g.resolve_record(gate['ctest'][key],root,stage)
checked=v.copy_qualification_test_inputs(stage/'generated/sirius/native_build_gate.json',stage);assert len(checked)==20
receipt={'kind':'reconstructed-exact-native-build-for-bounded-controls','source_revision':revision,'source_tree_clean':True,'original_ci_run':37922817020,'native_runtime_qualification_claimed':False,'ctest_registration_reconstructed':False,'source_root_checkout_eol_differences':source_eol,'whole_reconstructed_source_byte_gate_claimed':False,'canonical_and_consumed_inputs_checked':len(checked),'executables':{n:{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for n,p in executables.items()},'source_build_gate_sha256':hashlib.sha256((bundle/'native_build_gate.json').read_bytes()).hexdigest()}
(root/'out/native-wide-selection-review/reconstruction.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
