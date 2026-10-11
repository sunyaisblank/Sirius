"""Pin official decoder 0.1.6 headers/package and inspect its ELF metadata only."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import base64, hashlib, json, subprocess, tarfile, urllib.request
ROOT = Path.cwd(); W = Path(__file__).resolve().parent
REV = 'df604c2e6bc46bef23932d984b0b19b7d859bbe4'
def seal(p):
    data = p.read_bytes(); return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
def write(name, data):
    p = W/name; assert not p.exists(); p.write_bytes(data); return {'path': str(p.relative_to(ROOT)), **seal(p)}
def api(path):
    return subprocess.check_output(['gh', 'api', 'repos/ROCm/rocprof-trace-decoder/'+path], timeout=30)
ref_bytes = api('git/ref/tags/0.1.6'); ref = json.loads(ref_bytes)
assert ref['object'] == {'sha': REV, 'type': 'commit', 'url': 'https://api.github.com/repos/ROCm/rocprof-trace-decoder/git/commits/'+REV}
ref_record = write('decoder016-tag.json', ref_bytes)
release_bytes = api('releases/tags/0.1.6'); release = json.loads(release_bytes)
assert release['tag_name'] == '0.1.6' and not release['draft'] and not release['prerelease']
release_record = write('decoder016-release.json', release_bytes)
files = []
for source, name in [('include/rocprof_trace_decoder.h','decoder016-api.h'), ('include/trace_decoder_types.h','decoder016-types.h'), ('LICENSE','decoder016-LICENSE'), ('README.md','decoder016-README.md')]:
    meta = json.loads(api('contents/'+source+'?ref='+REV))
    assert meta['type']=='file' and meta['encoding']=='base64'
    data = base64.b64decode(meta['content']); assert len(data)==meta['size']
    blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest(); assert blob==meta['sha']
    files.append({**write(name,data),'repository':'ROCm/rocprof-trace-decoder','revision':REV,'repository_path':source,'git_blob':blob})
asset = next(a for a in release['assets'] if a['name']=='rocprof-trace-decoder-manylinux-2.28-0.1.6-Linux.tar.gz')
assert asset['size']==88697 and asset['digest']=='sha256:57a49c74db7dcd29433888ef5d408ca8e22223f2680c566855e4c3c68da7b36d'
assert asset['browser_download_url'].startswith('https://github.com/ROCm/rocprof-trace-decoder/releases/download/0.1.6/')
with urllib.request.urlopen(asset['browser_download_url'], timeout=30) as response:
    data = response.read(2*1024**2)
assert len(data)==asset['size'] and hashlib.sha256(data).hexdigest()==asset['digest'][7:]
package = write(asset['name'],data); members = []
lib_dir = ROOT/'bin/linux-gcc/qualification-tools/rocprof-trace-decoder-0.1.6'; assert not lib_dir.exists()
with tarfile.open(W/asset['name'],'r:gz') as archive:
    for member in archive.getmembers():
        assert member.isdir() or member.isfile()
        assert not Path(member.name).is_absolute() and '..' not in Path(member.name).parts
        members.append({'path':member.name,'bytes':member.size,'directory':member.isdir()})
    libs = [m for m in archive.getmembers() if m.isfile() and m.name.endswith('/librocprof-trace-decoder.so')]
    assert len(libs)==1 and libs[0].size<1024**2
    data = archive.extractfile(libs[0]).read(); assert len(data)==libs[0].size
lib_dir.mkdir(parents=True); lib = lib_dir/'librocprof-trace-decoder.so'; lib.write_bytes(data)
dynamic = write('decoder016-elf-dynamic.txt',subprocess.check_output(['readelf','-d','-V',str(lib)],timeout=10))
exports = write('decoder016-elf-exports.txt',subprocess.check_output(['readelf','--dyn-syms','--wide',str(lib)],timeout=10))
result = {'decoder_revision':REV,'release':release_record,'tag':ref_record,'files':files,'release_asset':{**package,'url':asset['browser_download_url'],'official_release_digest':asset['digest']},'members':members,'extracted_library':{'path':str(lib.relative_to(ROOT)),'archive_member':libs[0].name,**seal(lib)},'elf_dynamic':dynamic,'elf_exports':exports,'host_glibc':subprocess.check_output(['getconf','GNU_LIBC_VERSION'],text=True).strip(),'input_seals':{str(Path(__file__).relative_to(ROOT)):seal(Path(__file__))},'scope':'Verified official package, pinned public headers/legal terms, member-only extraction and ELF dynamic metadata inspection. No decoder execution, system install, binary disassembly or GPU work.'}
write('decoder016-preparation.json',(json.dumps(result,indent=2)+'\n').encode())
print(json.dumps({'revision':REV,'library':result['extracted_library'],'files':len(files),'package':package}))
