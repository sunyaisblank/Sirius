"""Retrieve the finite public hardware-name binding for observed counter groups."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess,urllib.request
ROOT=Path.cwd();W=Path(__file__).resolve().parent
revision='5b7a7544eac3ebc7f0f62ff878c69f8f893f3583';repo='GPUOpen-Tools/gpu_performance_api'
tree=json.loads(subprocess.check_output(['gh','api','repos/'+repo+'/git/trees/'+revision+'?recursive=1']))
assert not tree['truncated'];entries={v['path']:v for v in tree['tree']};records=[]
for block in ('cpf','sqwgp','ta','tcp','gl1c','gl2c','gcea'):
 origin='source/auto_generated/gpu_perf_api_counter_generator/gpa_hw_counter_gfx11_'+block+'.cc'
 obj=entries[origin];assert obj['type']=='blob' and obj['size']<131072
 url='https://raw.githubusercontent.com/'+repo+'/'+revision+'/'+origin
 with urllib.request.urlopen(url,timeout=30) as r:data=r.read(131073)
 assert len(data)==obj['size'] and hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==obj['sha']
 p=W/('gpa-hardware-'+block+'.cc');assert not p.exists();p.write_bytes(data)
 records.append({'path':str(p.relative_to(ROOT)),'repository_path':origin,'url':url,'git_blob':obj['sha'],'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
(W/'gpa-hardware-preparation.json').write_text(json.dumps({'revision':revision,'files':records},indent=2)+'\n');print(json.dumps({'files':len(records),'revision':revision}))
