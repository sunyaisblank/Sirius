"""Build sealed host inputs once; retain every other archive and link input."""
from pathlib import Path
import hashlib,json,re,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
doc=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
p=doc(W/'prepared-build.json');base=doc(R/p['base_P_result'])
def verify():
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==p['revision']
 assert not subprocess.check_output(['git','status','--porcelain']).strip()
 for name,expected in p['source_inputs'].items():assert pin(R/name)==expected,name
verify();assert not (W/'host-build-result.json').exists()
r={'passed':False,'commands':[],'source_inputs':p['source_inputs'],'source_proof':p['source_proof'],'scope':p['scope'],'new_cpp_compilations':0,'observer_compilations':0,'new_links':0,'new_shader_ISPC_provider_compilations':0,'same_actual_executable_for_both_provider_arms':True,'all3modules_once_in_both_executables':False,
 'include_path_relocation':base['include_path_relocation'],'spacing_probe':base['spacing_probe'],'unchanged_P_compute_source_proof':base['source_proof']}
def save():(W/'host-build-result.json').write_text(json.dumps(r,indent=2)+'\n')
def run(i):
 c=p['commands'][i];code=subprocess.run(c['argv'],cwd=c['cwd']).returncode;r['commands'].append({**c,'exit':code})
 if i==0:r['new_cpp_compilations']=1
 if i in (2,3):r['new_links']+=1
 save();assert code==0,(i,code)
run(0)
obj=R/p['object'];r['objects']=[{'argv':p['commands'][0]['argv'],'exit':0,'object':{'path':p['object'],**pin(obj)}}]
old=R/base['backend_archive'];archive=R/p['backend_archive'];assert not archive.exists();shutil.copyfile(old,archive);run(1)
members=subprocess.check_output(['/usr/bin/ar','t',str(old)],text=True).splitlines();assert members==subprocess.check_output(['/usr/bin/ar','t',str(archive)],text=True).splitlines()
different=[];unchanged={}
for name in members:
 a=subprocess.check_output(['/usr/bin/ar','p',str(old),name]);b=subprocess.check_output(['/usr/bin/ar','p',str(archive),name])
 if a!=b:different.append(name)
 else:unchanged[name]={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
assert different==['retained_integrator.cpp.o'];r.update(changed_archive_members=different,unchanged_archive_members=unchanged,backend_archive=p['backend_archive'],backend_archive_identity=pin(archive));save();run(2);run(3)
r['runtime_libraries']={};r['module_joins']={}
for key in ('executable','frame_executable'):
 path=R/p[key];r[key]=p[key];r[key+'_identity']=pin(path)
 text=subprocess.check_output(['ldd','-r',str(path)],text=True,stderr=subprocess.STDOUT);assert 'not found' not in text and 'undefined symbol' not in text
 for name in re.findall(r'(/[^\s()]+)',text):
  q=Path(name)
  if q.is_file():q=q.resolve();r['runtime_libraries'][str(q)]=pin(q)
for stage,v in base['module_joins'].items():
 q=R/v['module'];code=q.read_bytes();assert r['executable_identity'] and (R/p['executable']).read_bytes().count(code)==(R/p['frame_executable']).read_bytes().count(code)==1
 r['module_joins'][stage]={**v,'complete_occurrences':1,'executable_complete_occurrences':1,'frame_executable_complete_occurrences':1}
r.update(all3modules_once_in_both_executables=True,unchanged_other_link_inputs=base['unchanged_other_link_inputs'],frame_resources=base['frame_resources'],frame_resource_payload_bytes_copied=0,source_inputs_unchanged=True,original_backend_archive_unchanged=pin(old)==base['backend_archive_identity'],restored_private_headers=base['restored_private_headers'],restored_observer=base['restored_observer'])
verify();r['passed']=True;save();print(json.dumps({'passed':True,'new_cpp_compilations':1,'new_links':2,'changed_archive_members':different,'all3modules_once_in_both_executables':True}))
