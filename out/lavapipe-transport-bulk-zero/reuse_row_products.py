"""Restore only authenticated existing complete-row products after actual code admission."""
from pathlib import Path
import hashlib,json,os,subprocess,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent
def doc(p):return json.loads(p.read_text())
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
code=doc(W/'ordinary/ordinary-code-review.json');assert code['rows']['Transport']['remaining_optimized_IR_equivalent']
D=R/'attestations/qualification/bfa095ba/native-cpu-row-ownership-finite-outcome'
M=D/'MANIFEST.json';commit='9e01dfdca7798c79927827334798c72c2774eebc'
assert subprocess.check_output(['git','show',commit+':'+str(M.relative_to(R))])==M.read_bytes()
m=doc(M);old=R/'out/lavapipe-intergroup-row-owner/reused-row-products.json'
assert pin(old)=={'bytes':3212,'sha256':'a333b19b302b6de86a7069f995ef4d05eeb7d00a526485d390c89cbfd0d448c1'}
products=doc(old)['restored_actual_products'];archives={};restored=[]
for k,expected in products.items():
 v=m['actual_product_paths'][k];assert all(v[j]==expected[j] for j in ['bytes','sha256','mode'])
 p=R/k
 if not p.exists():
  blob=m['content_addressed_blobs'][v['sha256']];assert not blob['reused']
  a=D/blob['archive']
  if a not in archives:
   assert pin(a)==m['archives'][a.name]
   assert subprocess.check_output(['git','show',commit+':'+str(a.relative_to(R))])==a.read_bytes()
   archives[a]=pin(a)
  with zipfile.ZipFile(a) as z:b=z.read(blob['member'])
  assert len(b)==v['bytes'] and hashlib.sha256(b).hexdigest()==v['sha256']
  p.parent.mkdir(parents=True,exist_ok=True)
  with p.open('xb') as f:f.write(b)
  os.chmod(p,v['mode']);restored.append(k)
 assert pin(p)=={j:v[j] for j in ['bytes','sha256']} and p.stat().st_mode&0o777==v['mode']
matched=doc(R/'out/native-cpu-row-ownership/matched-result.json');exe=R/matched['executable']
assert pin(exe)['sha256']==matched['executable_sha256']=='7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'
compiler=doc(R/'out/native-cpu-row-ownership/compiler-result.json')
for kind,record in compiler['modules'].items():
 p=R/record['module'];assert pin(p)=={k:record[k] for k in ['bytes','sha256']}
 assert exe.read_bytes().count(p.read_bytes())==1,kind
result={'preservation_commit':commit,'source':'bfa095baf0ff9592cf648cfd49a73d5f6814c05c',
 'prior_manifest':{'path':str(M.relative_to(R)),**pin(M)},'restored_actual_products':products,
 'newly_restored_paths':restored,'authenticated_archives':{str(k.relative_to(R)):v for k,v in archives.items()},
 'three_complete_module_spans_present_once_in_final_ELF':True,'new_CPP_or_shader_compilations':0,
 'prior_results_reused_as_current_numerical_or_benefit_admission':False}
p=W/'reused-row-products.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
for name in ['compiler-result.json','matched-result.json']:
 p=W/name;assert not p.exists();original=R/'out/native-cpu-row-ownership'/name;record=doc(original)
 record['historical_result']={'path':str(original.relative_to(R)),**pin(original)}
 record['historical_commands']=record.get('commands',[]);record['commands']=[]
 record['compilation_performed_this_candidate']=False
 if name=='compiler-result.json':
  for module in record['modules'].values():module['compiled_this_candidate']=False
 else:
  record['new_cpp_compilations']=0;record['observer_compilations']=0
 record['reuse_scope']='Authenticated existing complete-row products restored exactly; original source/commands remain historical, zero new C++/shader/observer compilation and no old numerical/benefit admission.'
 p.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'actual_existing_products_verified':len(products),'paths_restored':restored,'new_compilations':0}))
