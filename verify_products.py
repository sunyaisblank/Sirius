from pathlib import Path
import json,hashlib,re,struct,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
source=json.loads((WORK/'source.json').read_text());assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==source['source_revision'];assert not subprocess.check_output(['git','status','--porcelain'])
old=json.loads((WORK/'baseline-products.json').read_text());baseline=json.loads((WORK/'baseline-arrays.json').read_text())['arrays']
def seal(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
for p in old['modules']:
 assert seal(ROOT/p['path'])=={k:p[k] for k in ['bytes','sha256']},p['path']
header=ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h';arrays={}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
 words=[int(q[:-1],0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)];assert len(words)==int(count) and name not in arrays
 data=struct.pack('<'+'I'*len(words),*words);arrays[name]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
assert len(arrays)==42 and len(baseline)==40
for name,identity in baseline.items():assert arrays[name]==identity,name
extra={'kTransportNarrowFmaShader':('retained_transport_narrow_fma.spv',87172,'f2755145d65fbc77cac634370a3a65143dada8d1665a6f1e3a99c689c4d00970'),'kEndpointNarrowFmaShader':('retained_endpoint_narrow_fma.spv',101108,'62b275566fa06b27872b7046da6ef721c383b2bb49087d6018965052241af611')}
assert set(arrays)-set(baseline)==set(extra)
for name,(filename,count,sha) in extra.items():
 assert arrays[name]=={'bytes':count,'sha256':sha}==seal(header.parent/filename)
 assert f'inline constexpr bool {name.replace("Shader","Available")} = true;' in header.read_text()
 query=subprocess.check_output(['ninja','-C','bin/linux-gcc','-t','query',str((header.parent/filename).relative_to(ROOT/'bin/linux-gcc'))],text=True);assert 'CUSTOM_COMMAND' in query
 (WORK/(filename+'.ninja-query.txt')).write_text(query)
assert len(list(header.parent.glob('*.spv')))==35
out={'source_revision':source['source_revision'],'baseline_revision':source['baseline_revision'],'whole_header_sha256':seal(header)['sha256'],'arrays':arrays,'original40_arrays_unchanged':True,'original33_modules_unchanged':True,'new2_match_independently_reviewed_compiler_products':True,'programs_unchanged':7,'optional_available':True,'cmake_byproducts_verified':2,'scope':'Integrated source/build/module identity only; no native execution, speed, frame or scientific/release qualification.'}
(WORK/'build-and-arrays.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['source_revision','original40_arrays_unchanged','original33_modules_unchanged','new2_match_independently_reviewed_compiler_products','whole_header_sha256']}))
