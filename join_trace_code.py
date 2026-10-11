"""Join saved driver code to preserved original source-bound query products."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,re,struct,subprocess
ROOT=Path.cwd();W=Path(__file__).resolve().parent
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def identity(b):return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def git(*a):return subprocess.check_output(['git',*a]).decode().strip()
head=git('rev-parse','HEAD');assert head==git('rev-parse','@{u}') and not git('status','--porcelain')
assert git('rev-parse','HEAD^{tree}')==git('rev-parse','b81e061^{tree}')
plan=json.loads((W/'counter-view-plan.json').read_text());trace=ROOT/plan['protected_trace']['path'];assert seal(trace)=={k:plan['protected_trace'][k] for k in ('bytes','sha256')}
index_blob=subprocess.check_output(['git','show','88d3a0f2f001bbf1037902e0bf33fd9111c5f159:spm-index-result.json']);index=json.loads(index_blob)
index_path=W/'code-join-index.json'
if index_path.exists():assert index_path.read_bytes()==index_blob
else:index_path.write_bytes(index_blob)
q=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code';prep=q/'preparation-query.json';receipt=json.loads(prep.read_text());old=receipt['revision']
original_elf=q/'query-payloads/transport.info1';endpoint_isa=q/'endpoint-payloads/endpoint.info2'
header=ROOT/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h';assert seal(header)['sha256']==receipt['source_bindings'][str(header.relative_to(ROOT))]
backend=ROOT/'bin/windows-msvc/native-b81e061-executed-cost/tests/backend/Release/sirius_backend_tests.exe';pe=backend.read_bytes()
paths=[Path(__file__),W/'counter-view-plan.json',trace,index_path,prep,original_elf,endpoint_isa,header,backend]
modules=[]
for m in receipt['modules']:
 p=q/m['path'];b=p.read_bytes();assert identity(b)=={k:m[k] for k in ('bytes','sha256')}
 offset=pe.find(b);assert offset>=0 and pe.find(b,offset+1)==-1
 modules.append({'stage':m['stage'],**identity(b),'current_native_backend_offset':offset,'complete_original_module_bytes_equal':True});paths.append(p)
shared=[ROOT/name for name in ('src/sirius/kernels/retained_triple.slang','src/sirius/kernels/retained_pair.slang','src/sirius/kernels/retained_scalar.slang','src/sirius/kernels/retained_scratch.slang','src/sirius/kernels/portable_binary32.h','src/sirius/kernels/retained_transport.slang','src/sirius/kernels/retained_endpoint.slang','scripts/build-retained-kernels.py','src/sirius/kernels/retained_program.py')]
kernel_sources=[]
for p in shared:
 name=str(p.relative_to(ROOT))
 prior=subprocess.check_output(['git','show',old+':'+name],stderr=subprocess.DEVNULL)
 assert p.read_bytes()==prior,name;kernel_sources.append({'path':name,**seal(p),'exact_original_query_revision':old});paths.append(p)
def sections(b):
 assert b[:6]==b'\x7fELF\x02\x01'
 offset=struct.unpack_from('<Q',b,40)[0];size,n,names=struct.unpack_from('<3H',b,58);assert size==64 and 0<n<=256 and names<n and offset+n*size<=len(b)
 rows=[struct.unpack_from('<IIQQQQIIQQ',b,offset+i*size) for i in range(n)];r=rows[names];strings=b[r[4]:r[4]+r[5]];out={}
 for i,r in enumerate(rows):
  assert r[0]<len(strings) and r[4]+r[5]<=len(b)
  name=strings[r[0]:].split(b'\0',1)[0].decode();assert name not in out
  out[name]={'index':i,'type':r[1],'flags':r[2],'link':r[6],'offset':r[4],'data':b[r[4]:r[4]+r[5]]}
 return out
b=trace.read_bytes();objects={}
for c in index['chunks']:
 if c['name']!='CodeObject':continue
 data=b[c['data_offset']:c['data_offset']+c['data_bytes']];assert identity(data)=={'bytes':c['data_bytes'],'sha256':c['data_sha256']}
 pci=struct.unpack_from('<I',b,c['header_offset'])[0];lo,hi=struct.unpack_from('<2Q',b,c['header_offset']+8)
 objects[c['ordinal']]={'chunk':c,'data':data,'sections':sections(data),'pci_id':pci,'internal_lower':f'{lo:016x}','internal_upper':f'{hi:016x}'}
t=objects[4];e=objects[3];original=sections(original_elf.read_bytes());ts=t['sections'];es=e['sections']
assert {k for k,v in ts.items() if v['flags']&2}=={'.text','.rodata'}=={k for k,v in original.items() if v['flags']&2}
joined={}
for name in ('.text','.rodata','.rel.text'):
 assert ts[name]['data']==original[name]['data'];joined[name]=identity(ts[name]['data'])
def symbols(s):
 table=s['.symtab'];assert table['type']==2 and len(table['data'])%24==0
 names={v['index']:k for k,v in s.items()};strings=s[names[table['link']]]['data'];out=[]
 for off in range(0,len(table['data']),24):
  name,info,other,section,value,size=struct.unpack_from('<IBBHQQ',table['data'],off)
  assert name<len(strings)
  out.append({'name':strings[name:].split(b'\0',1)[0].decode(),'info':info,'other':other,'section':names[section] if section in names else section,'value':value,'size':size})
 return out
trace_symbols=symbols(ts);original_symbols=symbols(original)
assert len(trace_symbols)==9 and len(original_symbols)==10
assert trace_symbols==original_symbols[:len(trace_symbols)]
assert original_symbols[9]['name']=='llvmir' and original_symbols[9]['section']=='.AMDGPU.comment.llvmir'
def normalize_note(note):
 result=bytearray(note)
 for key in (b'.internal_pipeline_hash',b'.128_bit_cache_hash'):
  assert note.count(key)==1;offset=note.index(key)+len(key);assert note[offset]==0x92;offset+=1
  for _ in range(2):
   assert note[offset]==0xcf;result[offset+1:offset+9]=bytes(8);offset+=9
 return bytes(result)
assert normalize_note(ts['.note']['data'])==normalize_note(original['.note']['data'])
isa=endpoint_isa.read_bytes();assert hashlib.sha256(isa).hexdigest()=='2b706dc8e57107210a1422583627764d0fac1413f63d83416d15f6ea7a2834ed'
pattern=rb'^\s*([a-z][a-z0-9_]*)(?=\s).*?;\s+([a-f0-9]+(?:\s+[a-f0-9]+)*)\s*$'
instructions=list(re.finditer(pattern,isa,re.M));encoded=b''.join(struct.pack('<I',int(word,16)) for m in instructions for word in m[2].split())
assert len(instructions)==130811 and len(encoded)==667368 and es['.text']['data'][:len(encoded)]==encoded
sym=es['.symtab']['data'];strings=es['.strtab']['data'];main=[]
for off in range(0,len(sym),24):
 name,info,other,section,value,size=struct.unpack_from('<IBBHQQ',sym,off)
 if strings[name:].split(b'\0',1)[0]==b'_amdgpu_cs_main':main.append((value,size))
assert main==[(0,len(encoded))]
tail=es['.text']['data'][len(encoded):];assert tail==struct.pack('<I',0xbf9f0000)*102
correlation=next(c for c in index['chunks'] if c['name']=='PsoCorrelation');assert correlation['version']==3 and correlation['data_bytes']==288
correlations=[]
for off in range(correlation['data_offset'],correlation['data_offset']+288,96):
 pci=struct.unpack_from('<I',b,off)[0];api,lo,hi=struct.unpack_from('<3Q',b,off+8)
 correlations.append({'pci_id':pci,'api_pso':f'{api:016x}','internal_lower':f'{lo:016x}','internal_upper':f'{hi:016x}'})
for o,api in ((t,'b4bad5049f1d989e'),(e,'a08a2f262eca2b69')):
 assert {'pci_id':o['pci_id'],'api_pso':api,'internal_lower':o['internal_lower'],'internal_upper':o['internal_upper']} in correlations
inputs={str(p.relative_to(ROOT)):seal(p) for p in paths}
result={'source_revision':head,'input_seals':inputs,'source_bound_query_revision':old,'current_original_modules':modules,'kernel_source_joins':kernel_sources,'Transport':{'code_object':{k:t['chunk'][k] for k in ('ordinal','data_offset','data_bytes','data_sha256')},'allocated_sections_and_relocations_exact':joined,'symbols_semantically_exact':trace_symbols,'raw_symbol_table_equality':False,'extra_baseline_nonallocated_symbol':original_symbols[9],'pal_note_equal_after_only_internal_and_cache_hash_normalization':True,'internal_pipeline_lower':t['internal_lower'],'internal_pipeline_upper':t['internal_upper']},'Endpoint':{'code_object':{k:e['chunk'][k] for k in ('ordinal','data_offset','data_bytes','data_sha256')},'complete_declared_entry_encoding_exact':identity(encoded),'instructions':len(instructions),'uncompared_text_tail':identity(tail),'tail_words':'102 x bf9f0000','scope':'Entire declared entry extent equals source-bound printed ISA. No archived whole Endpoint ELF comparator; unprinted tail/other sections whole-byte equality unproved.'},'pso_correlations':correlations,'scope':'Exact saved code/source-product joins, with stated Endpoint limitation. Actual selected-event API/internal hashes still required. No raw-SPM timestamp conversion, exclusive counter ownership, optimization benefit, rendered frame or full qualification.'}
for name,value in inputs.items():assert seal(ROOT/name)==value
p=W/'trace-code-join.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'Transport_allocated_code_constants_relocations_exact':True,'Endpoint_complete_entry_bytes_exact':len(encoded),'current_original_modules':len(modules),'current_kernel_source_joins':len(kernel_sources)}))
