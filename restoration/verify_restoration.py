"""Verify exact original source/products after the rejected eight-file trial."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,re,struct,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent;TRIAL=WORK.parent
def read(p):return json.loads(p.read_text())
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
source=read(WORK/'source.json');head=source['source_revision'];baseline=source['baseline_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'])
assert subprocess.check_output(['git','rev-parse',head+'^{tree}'])==subprocess.check_output(['git','rev-parse',baseline+'^{tree}'])
old=read(TRIAL/'baseline-products.json');arrays0=read(TRIAL/'baseline-arrays.json')['arrays']
assert old['revision']==baseline and len(old['modules'])==33 and len(arrays0)==40
for v in old['modules']:
    assert seal(ROOT/v['path'])=={k:v[k] for k in ('bytes','sha256')},v['path']
header=ROOT/old['header']['path']
assert seal(header)=={k:old['header'][k] for k in ('bytes','sha256')}
assert header.read_bytes()==(TRIAL/'baseline-retained_kernels.h').read_bytes()
arrays={}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
    words=[int(q[:-1],0) for q in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)]
    assert len(words)==int(count) and name not in arrays
    b=struct.pack('<'+'I'*len(words),*words);arrays[name]={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
assert arrays==arrays0 and len(arrays)==40
declared={Path(v['path']).name for v in old['modules']}
obsolete={'retained_transport_narrow_fma.spv','retained_endpoint_narrow_fma.spv'}
physical={p.name for p in header.parent.glob('*.spv')}
assert physical==declared|obsolete
for name in obsolete:
    q=subprocess.run(['ninja','-C','bin/linux-gcc','-t','query',str((header.parent/name).relative_to(ROOT/'bin/linux-gcc'))],capture_output=True,text=True)
    assert q.returncode!=0 and 'unknown target' in q.stderr,(name,q.stdout,q.stderr)
record={'source_revision':head,'baseline_revision':baseline,'whole_header_sha256':seal(header)['sha256'],
    'arrays':arrays,'tracked_tree_exact_baseline':True,'original40_arrays_unchanged':True,
    'original33_modules_unchanged':True,'whole_header_exact_baseline':True,'programs_unchanged':7,
    'obsolete_nondeclared_byproducts_held_pending_cleanup':sorted(obsolete),
    'scope':'Exact restored source/build/module identity only. Original two obsolete trial SPVs remain protected until final cleanup. No fresh native/scientific/render/release qualification.'}
(WORK/'build-and-arrays.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ('source_revision','tracked_tree_exact_baseline','whole_header_exact_baseline','original40_arrays_unchanged','original33_modules_unchanged')}))
