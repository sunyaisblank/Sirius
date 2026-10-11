"""Pin the public code-load record schema and its actual producer."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess,urllib.request
ROOT=Path.cwd();W=Path(__file__).resolve().parent
repo='llvm/llvm-project';revision='3b5b5c1ec4a3095ab096dd780e84d7ab81f3d7ff';records=[]
for origin,name in [('llvm/lib/Target/AMDGPU/MCTargetDesc/AMDGPUInstPrinter.cpp','llvm-instruction-printer.cpp'),('llvm/lib/Target/AMDGPU/Disassembler/AMDGPUDisassembler.cpp','llvm-disassembler.cpp')]:
 obj=json.loads(subprocess.check_output(['gh','api','repos/'+repo+'/contents/'+origin+'?ref='+revision]));assert obj['type']=='file' and obj['size']<1048576
 url='https://raw.githubusercontent.com/'+repo+'/'+revision+'/'+origin
 with urllib.request.urlopen(url,timeout=30) as r:data=r.read(1048577)
 assert len(data)==obj['size'] and hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==obj['sha']
 p=W/name;assert not p.exists();p.write_bytes(data);records.append({'path':str(p.relative_to(ROOT)),'url':url,'revision':revision,'git_blob':obj['sha'],'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
(W/'llvm-branch-preparation.json').write_text(json.dumps({'files':records},indent=2)+'\n');print(json.dumps({'pinned_files':len(records)}))
