"""Pin primary consumer sources and extract only the verified offline decoder."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess,tarfile,urllib.request
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(data):return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def gh(path):return json.loads(subprocess.check_output(['gh','api',path]))
def fetch(url):
 with urllib.request.urlopen(url,timeout=30) as r:data=r.read(131073)
 assert len(data)<=131072
 return data
sources=[]
for repo,origin,name in [('ROCm/rocprofiler-sdk','source/lib/rocprofiler-sdk/thread_trace/decode.cpp','sdk-decode.cpp'),('ROCm/aqlprofile','src/core/threadtrace.cpp','aql-threadtrace.cpp')]:
 revision=gh('repos/'+repo+'/commits/rocm-7.1.0')['sha']
 obj=gh('repos/'+repo+'/contents/'+origin+'?ref='+revision)
 data=fetch('https://raw.githubusercontent.com/'+repo+'/'+revision+'/'+origin)
 assert obj['type']=='file' and len(data)==obj['size']
 assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==obj['sha']
 p=W/name;assert not p.exists();p.write_bytes(data)
 sources.append({'path':str(p.relative_to(ROOT)),'repository':repo,'revision':revision,'repository_path':origin,'git_blob':obj['sha'],**seal(data)})
prep=json.loads((W/'decoder-preparation.json').read_text());asset=prep['release_asset']
p=ROOT/asset['path'];assert seal(p.read_bytes())=={k:asset[k] for k in ('bytes','sha256')}
regular=[v for v in prep['members'] if not v['directory']];assert len(regular)==1
member=regular[0];assert member['bytes']==199792
with tarfile.open(p,'r:gz') as t:
 m=t.getmember(member['path']);assert m.isfile() and not m.issym() and not m.islnk()
 data=t.extractfile(m).read();assert len(data)==member['bytes']
target=ROOT/'bin/linux-gcc/qualification-tools/rocprof-trace-decoder-0.1.3/librocprof-trace-decoder.so'
assert subprocess.check_output(['git','check-ignore',str(target.relative_to(ROOT))],text=True).strip()==str(target.relative_to(ROOT))
assert not target.exists();target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
dependency=subprocess.run(['readelf','-d',str(target)],capture_output=True,text=True,check=True)
(W/'decoder-elf-dynamic.txt').write_text(dependency.stdout)
result={'input_seals':{str(Path(__file__).relative_to(ROOT)):seal(Path(__file__).read_bytes()),str(p.relative_to(ROOT)):seal(p.read_bytes())},'primary_sources':sources,'extracted_library':{'path':str(target.relative_to(ROOT)),'archive_member':member['path'],**seal(data)},'host_glibc':subprocess.check_output(['getconf','GNU_LIBC_VERSION'],text=True).strip(),'elf_dynamic':{'path':str((W/'decoder-elf-dynamic.txt').relative_to(ROOT)),**seal(dependency.stdout.encode())},'scope':'Verified member-only extraction and static ELF inspection under existing ignored preset. No library execution, system install, HSA/ROCm/GPU collection or product change.'}
(W/'cpu-decoder-preparation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
