"""Complete producer/module frame conditions before integrated execution."""
from pathlib import Path
import re,struct,hashlib,json
r=Path.cwd();w=Path(__file__).resolve().parent;b=r/'bin/linux-gcc/cooperative-dense-inputs';base=r/'bin/linux-gcc/src/sirius/backend/retained';checked=r/'attestations/software-vulkan/dcd1368/sealed-stage-seven-owner/files/bin/linux-gcc/sealed-stage-seven-owner/products'
def arrays(p):
 d={}
 for size,n,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',p.read_text(),re.S):
  values=[int(v[:-1],0) for v in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)];assert len(values)==int(size);data=struct.pack('<%dI'%len(values),*values);d[n]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
 return d
def module(p):
 data=p.read_bytes();a=struct.unpack('<%dI'%(len(data)//4),data);assert a[0]==0x07230203;const={};joins=[];caps=[];types=[];modes=[];i=5
 while i<len(a):
  size,op=a[i]>>16,a[i]&0xffff;assert size and i+size<=len(a);args=a[i+1:i+size]
  if op==43 and size==4:const[args[1]]=args[2]
  if op==224:joins.append(args)
  if op==17:caps.append(args)
  if op==21:types.append(args[1:])
  if op==16:modes.append(args[1:])
  i+=size
 return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'capabilities':caps,'integer_types':sorted(set(types)),'barriers':sorted(tuple(const[v] for v in x) for x in joins),'execution_modes':modes}
old=arrays(base/'retained_kernels.h');new=arrays(b/'products/retained_kernels.h');assert old.keys()==new.keys() and len(old)==40;changed_arrays={'kDensePortableShader','kDensePortableFp64Shader'};assert all(old[n]==new[n] for n in old if n not in changed_arrays);assert {n for n in old if old[n]!=new[n]}==changed_arrays
for name,d in [('baseline-arrays.json',old),('candidate-arrays.json',new)]:(w/name).write_text(json.dumps({'arrays':d},indent=2)+'\n')
records=[]
for p in sorted((b/'products').glob('*.spv')):
 q=base/p.name
 if not q.exists():q=checked/p.name
 assert q.is_file(),q;a=module(q);c=module(p)
 for key in ('capabilities','integer_types','execution_modes'):assert a[key]==c[key],(p,key)
 if p.name in ('retained_dense_portable.spv','retained_dense_portable_fp64.spv'):
  from collections import Counter
  previous=Counter(map(tuple,a['barriers']));current=Counter(map(tuple,c['barriers']));assert previous<=current and current-previous==Counter({(2,2,264):2}),(p,'original barriers plus exactly two group joins',previous,current)
 else:assert a==c,(p,'unaffected module byte identity')
 records.append({'module':p.name,'baseline_path':str(q.relative_to(r)),'baseline':a,'candidate':c,'same_bytes':a['sha256']==c['sha256']})
assert len(records)==35
result={'passed':True,'complete_products':35,'all7program_arrays_and16native_modules_exact':True,'capabilities_integer_types_layouts_exact_original_barriers_preserved_plus2Dense_group_joins':True,'arrays':{'same':[n for n in old if old[n]==new[n]],'changed':[n for n in old if old[n]!=new[n]]},'records':records};(w/'generated-products-check.json').write_text(json.dumps(result,indent=2)+'\n');print({'products':35,'unchanged_modules':sum(x['same_bytes'] for x in records),'unchanged_arrays':len(result['arrays']['same'])})
