"""Pin the public code-load record schema and its actual producer."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess,urllib.request
ROOT=Path.cwd();W=Path(__file__).resolve().parent
repo='ROCm/llvm-project';revision='1b0eada6b0ee93e2e694c8c146d23fca90bc11c5';records=[]
for origin,name in [('llvm/LICENSE.TXT','rocm-llvm-LICENSE.txt')]:
 obj=json.loads(subprocess.check_output(['gh','api','repos/'+repo+'/contents/'+origin+'?ref='+revision]));assert obj['type']=='file' and obj['size']<131072
 url='https://raw.githubusercontent.com/'+repo+'/'+revision+'/'+origin
 with urllib.request.urlopen(url,timeout=30) as r:data=r.read(131073)
 assert len(data)==obj['size'] and hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==obj['sha']
 p=W/name;assert not p.exists();p.write_bytes(data);records.append({'path':str(p.relative_to(ROOT)),'url':url,'revision':revision,'git_blob':obj['sha'],'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
(W/'rocm-license-preparation.json').write_text(json.dumps({'files':records},indent=2)+'\n');print(json.dumps({'pinned_files':len(records)}))
