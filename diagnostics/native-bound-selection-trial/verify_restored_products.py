"""Verify restored generated payloads against the sealed accepted baseline."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,re,struct,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
ACCEPTED='2053ef79c90752d37b0075b8137cfb9b571b6820'
def read(p):return json.loads(p.read_text())
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
def git(*a):return subprocess.check_output(['git',*a],text=True,timeout=30).strip()
assert WORK==ROOT/'out/native-bound-selection-trial'
assert git('rev-parse','HEAD')==ACCEPTED==git('rev-parse','@{u}') and not git('status','--porcelain')
build=read(WORK/'restored-source-build.json');assert build['accepted'] and build['source_revision']==ACCEPTED and not build['errors']
for name,expected in build['input_seals'].items():assert seal(Path(name))==expected
baseline=read(WORK/'baseline-arrays.json');assert baseline['source_revision']==ACCEPTED
old=WORK/'baseline/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
header=ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert seal(old)=={'bytes':baseline['header_bytes'],'sha256':baseline['header_sha256']}
assert header.read_bytes()==old.read_bytes()
arrays={};payloads={}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
 words=[int(x[:-1],0) for x in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)]
 assert name not in arrays and len(words)==int(count)
 data=struct.pack('<'+'I'*len(words),*words);payloads[name]=data
 arrays[name]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
assert len(arrays)==40 and arrays==baseline['arrays']
modules={}
for stem,array in (('retained_transport','kTransportShader'),('retained_endpoint','kEndpointShader')):
 p=header.parent/(stem+'.spv');assert p.read_bytes()==payloads[array];modules[stem]=seal(p)
source={}
for name in ('src/sirius/kernels/retained_pair.slang','scripts/build-retained-kernels.py'):
 p=ROOT/name;assert p.read_bytes()==(WORK/'baseline'/name).read_bytes();source[name]=seal(p)
record={'accepted':True,'source_revision':ACCEPTED,'header':seal(header),'complete_arrays_matched':40,'arrays':arrays,'modules':modules,'source':source,'build_receipt':seal(WORK/'restored-source-build.json'),'reader':seal(Path(__file__)),'scope':'Restored source and actual local generated payloads equal accepted baseline; source/build only, no new native run or scientific/release qualification.'}
p=WORK/'restored-products.json';assert not p.exists();p.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('accepted','source_revision','header','complete_arrays_matched','modules','scope')}))
