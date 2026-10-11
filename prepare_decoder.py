"""Fetch pinned public headers and one official offline decoder package."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,subprocess,tarfile,urllib.request
ROOT=Path.cwd();W=Path(__file__).resolve().parent
revision='0af1dbcd285946a84df0df7d3118934becc76441'
repo='ROCm/rocprof-trace-decoder'
def gh(path):return json.loads(subprocess.check_output(['gh','api',path]))
def fetch(url,limit):
 with urllib.request.urlopen(url,timeout=30) as response:
  data=response.read(limit+1);assert len(data)<=limit
 return data
def seal(data):return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
tree=gh('repos/'+repo+'/git/trees/'+revision+'?recursive=1');entries={v['path']:v for v in tree['tree']}
files=[]
for origin,target in (('include/rocprof_trace_decoder.h','rocprof_trace_decoder.h'),('include/trace_decoder_types.h','trace_decoder_types.h'),('LICENSE','decoder-LICENSE'),('README.md','decoder-README.md')):
 e=entries[origin];assert e['type']=='blob' and e['size']<32768
 url='https://raw.githubusercontent.com/'+repo+'/'+revision+'/'+origin
 data=fetch(url,32768);assert len(data)==e['size'] and hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==e['sha']
 p=W/target;assert not p.exists();p.write_bytes(data);files.append({'path':str(p.relative_to(ROOT)),'url':url,'revision':revision,'repository_path':origin,'git_blob':e['sha'],**seal(data)})
gpa_revision='5b7a7544eac3ebc7f0f62ff878c69f8f893f3583';origin='source/codegen/public_counter_definitions_gfx11.txt'
url='https://raw.githubusercontent.com/GPUOpen-Tools/gpu_performance_api/'+gpa_revision+'/'+origin;data=fetch(url,65536)
assert len(data)==51492 and hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()=='a1d289bda4b021d99f66abc5717d5e7b2b7adaaa'
p=W/'gpa-public-counter-definitions-gfx11.txt';assert not p.exists();p.write_bytes(data);files.append({'path':str(p.relative_to(ROOT)),'url':url,'revision':gpa_revision,'repository_path':origin,'git_blob':'a1d289bda4b021d99f66abc5717d5e7b2b7adaaa',**seal(data)})
release=gh('repos/'+repo+'/releases/tags/0.1.3');assert release['tag_name']=='0.1.3'
(W/'decoder-release.json').write_text(json.dumps(release,indent=2)+'\n')
name='rocprof-trace-decoder-ubuntu-24.04-0.1.3-Linux.tar.gz';matches=[v for v in release['assets'] if v['name']==name];assert len(matches)==1
asset=matches[0];assert asset['size']==82987 and asset['digest']=='sha256:d2c0a4147a7ac7a179da89fd5d7370cbfadb0a67c7ced9cd3ec27109d9c2f3c8'
data=fetch(asset['browser_download_url'],1048576);assert len(data)==asset['size'] and hashlib.sha256(data).hexdigest()==asset['digest'][7:]
p=W/name;assert not p.exists();p.write_bytes(data)
with tarfile.open(p,'r:gz') as archive:
 members=[]
 for m in archive.getmembers():
  assert not m.issym() and not m.islnk() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts and (m.isdir() or m.isfile())
  members.append({'path':m.name,'bytes':m.size,'directory':m.isdir()})
result={'decoder_revision':revision,'files':files,'release_asset':{'path':str(p.relative_to(ROOT)),'url':asset['browser_download_url'],'official_release_digest':asset['digest'],**seal(data)},'members':members,'scope':'Pinned public header/license/counter-definition retrieval and verified official package inventory only. No package install, code/library execution, GPU/HSA/ROCm collection, driver/environment mutation or decoder success claim.'}
(W/'decoder-preparation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'files':len(files),'asset':result['release_asset'],'members':members}))
