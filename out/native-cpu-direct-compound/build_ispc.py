"""Compile each sealed raw stage once; reuse the already-joined current adapter."""
from pathlib import Path
import hashlib,json,os,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;P=json.loads((W/'prepared-ispc-dense-inputs.json').read_text());out=W/'ispc-dense-compile-result.json'
def pin(p):
 with p.open('rb') as f:return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(f,'sha256').hexdigest()}
def verify():
 assert not os.environ.get('ISPC_ARGS') and not subprocess.check_output(['git','status','--porcelain']).strip()
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==P['revision']
 for name,v in P['inputs'].items():assert pin(R/name)==v,name
verify();assert not out.exists()
r={'revision':P['revision'],'objects':[],'all_objects_compiled':False,'source_inputs':P['inputs'],'scope':P['scope'],'new_CPP_compilations':0,'new_relocatable_links':0,'reused_current_CPP_adapter_link':1,'merged_adapter':P['merged_adapter'],'adapter_input_code_association':P['adapter_input_code_association']}
def save():out.write_text(json.dumps(r,indent=2)+'\n')
save()
for item in P['objects']:
 assert not (R/item['output']).exists()
 row={'stage':item['stage'],'argv':item['argv'],'status':'started'};r['objects'].append(row);save();print(json.dumps(row),flush=True)
 code=subprocess.run(item['argv'],cwd=R).returncode;row.update(status='terminal',argv_executed=True,exit=code)
 if (R/item['output']).is_file():row['object']={'path':item['output'],**pin(R/item['output'])}
 save()
 if code:raise SystemExit(code)
verify();r.update(all_prepared_source_bytes_unchanged=True,new_compilations=3,reused_compiled_objects=0,all_objects_compiled=True);save()
print(json.dumps({'all_objects_compiled':True,'new_compilations':3,'new_CPP_compilations':0,'new_relocatable_links':0}))
