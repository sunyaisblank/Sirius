"""Verify restored accepted source and actual read-only production payloads."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, struct, subprocess
ROOT = Path.cwd();WORK = Path(__file__).resolve().parent
assert ROOT.resolve() == Path('/home/astra/.project/Sirius') and WORK == ROOT/'out/native-division-first-correction-screen'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
def read(p): return json.loads(p.read_text())
def seal(p):
    data = p.read_bytes();return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=20).strip()
review = read(WORK/'restored-product-reader-review.json')
assert review['accepted'] and review['reader_seal'] == seal(Path(__file__))
baseline = read(WORK/'baseline-source.json');build = read(WORK/'restored-source-build.json')
assert baseline['revision'] == build['source_revision'] == BASELINE and build['accepted']
assert git('rev-parse','HEAD') == BASELINE == git('rev-parse','@{u}') and not git('status','--porcelain')
assert set(git('ls-files').splitlines()) == set(baseline['tracked_source_seals'])
for name,expected in baseline['tracked_source_seals'].items():assert seal(ROOT/name) == expected,name
assert not (ROOT/'tests/backend/shader_clock_observer.h').exists()
header = ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
assert seal(header) == baseline['production_header']
paths = [ROOT/name for name in baseline['tracked_source_seals']]
paths.extend((Path(__file__),WORK/'restored-product-reader-review.json',WORK/'baseline-source.json',WORK/'restored-source-build.json',header))
executables = {name:ROOT/relative for name,relative in {
    'app':'bin/linux-gcc/src/sirius/app/sirius',
    'backend_tests':'bin/linux-gcc/tests/backend/sirius_backend_tests'}.items()}
paths.extend(executables.values())
transport = WORK/'runtime/candidate/retained_transport.spv'
endpoint = WORK/'runtime/candidate/retained_endpoint.spv'
paths.extend((transport,endpoint))
inputs = {str(path):seal(path) for path in paths}
arrays = {}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
    words = [int(word[:-1],0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)]
    assert len(words) == int(count) and name not in arrays
    arrays[name] = struct.pack('<'+'I'*len(words),*words)
assert len(arrays) == 40
products = {}
for name,path in executables.items():
    data = path.read_bytes();assert data[:6] == b'\x7fELF\x02\x01'
    table = struct.unpack_from('<Q',data,0x28)[0]
    stride,count = struct.unpack_from('<HH',data,0x3a)
    assert stride == 64 and count > 0 and table+stride*count <= len(data)
    sections = []
    for i in range(count):
        fields = struct.unpack_from('<IIQQQQIIQQ',data,table+i*stride)
        kind,flags,address,offset,length = fields[1:6]
        if kind == 1 and flags & 2 and not flags & 5:
            assert offset+length <= len(data)
            sections.append({'section_index':i,'file_offset':offset,'address':address,'bytes':length,'flags':flags})
    assert sections
    bindings = {}
    for array,payload in arrays.items():
        matches = []
        for section in sections:
            begin = section['file_offset'];pos = data[begin:begin+section['bytes']].find(payload)
            if pos >= 0: matches.append({**section,'payload_file_offset':begin+pos,'payload_address':section['address']+pos})
        assert matches,array
        bindings[array] = {'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest(),'readonly_ELF_bindings':matches}
    assert transport.read_bytes() not in data and endpoint.read_bytes() not in data
    products[name] = {'executable':seal(path),'production_payloads':bindings,'rejected_candidate_modules_absent':True}
assert {str(path):seal(path) for path in paths} == inputs
assert git('rev-parse','HEAD') == BASELINE == git('rev-parse','@{u}') and not git('status','--porcelain')
result = {'accepted':True,'source_revision':BASELINE,'input_seals':inputs,'source_paths_exact':len(baseline['tracked_source_seals']),
    'production_header':seal(header),'production_arrays':40,'products':products,
    'scope':'Restored exact accepted source and unchanged whole production payloads in read-only app/backend ELF sections, with rejected candidate modules absent. Stored payload presence/source-build only; full numerical/rendering/scientific/release qualification remains open.'}
path = WORK/'restored-products.json';assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('accepted','source_revision','source_paths_exact','production_header','production_arrays')}))
