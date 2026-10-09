#!/usr/bin/env python3
"""Offline identity/layout verification and syntactic call-site accounting."""
import sys
sys.dont_write_bytecode=True
import collections, functools, hashlib, importlib.util, json, re, struct, subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PRESET=ROOT/'bin/linux-gcc/src/sirius/backend/retained'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
baseline=json.loads((HERE/'baseline.json').read_text())
changed=[]
for name,expected in baseline['original_source_sha256'].items():
 assert sha(ROOT/name)==expected,('tracked baseline changed',name)
 if sha(HERE/name)!=expected:changed.append(name)
stems=['camera','ray_camera','transport','endpoint','dense','initialize']
expected_changes={'scripts/build-retained-kernels.py','src/sirius/kernels/retained_triple.slang','src/sirius/kernels/retained_pair.slang',*[f'src/sirius/kernels/retained_{stem}.slang' for stem in stems]}
assert set(changed)==expected_changes,changed
for name in expected_changes-{'scripts/build-retained-kernels.py'}:
 candidate=(HERE/name).read_text();original=(ROOT/name).read_text()
 pattern=r'#ifdef SIRIUS_RETAINED_PORTABLE\n.*?#else\n(.*?)#endif\n'
 native,n=re.subn(pattern,lambda m:m[1],candidate,flags=re.S)
 assert n==(4 if name.endswith('retained_triple.slang') else 1),(name,n)
 native,nbar=re.subn(r'#ifndef SIRIUS_RETAINED_PORTABLE\n(.*?)#endif\n',lambda m:m[1],native,flags=re.S)
 native=native.replace('SIRIUS_RETAINED_EXECUTION_LANES','SIRIUS_RETAINED_LANES')
 assert native==original,('native branch text differs',name)
assert sha(HERE/'src/sirius/kernels/retained_scratch.slang')==sha(ROOT/'src/sirius/kernels/retained_scratch.slang')
spec=importlib.util.spec_from_file_location('retained_program',HERE/'src/sirius/kernels/retained_program.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
headers=[(p.read_text()) for p in [HERE/'modules/retained_kernels.h',PRESET/'retained_kernels.h']]
def words(text,name):
 match=re.search(r'inline constexpr std::array<std::uint32_t, (\d+)> '+name+r'\{\{(.*?)\}\};',text,re.S)
 assert match,name
 result=[int(v) for v in re.findall(r'(\d+)u',match[2])]
 assert len(result)==int(match[1]);return result

def call_counts(path):
 data=path.read_bytes();ws=struct.unpack('<'+str(len(data)//4)+'I',data)
 names={};functions={};calls=collections.Counter();entry=None;current=None;i=5
 while i<len(ws):
  count,opcode=ws[i]>>16,ws[i]&65535;assert count and i+count<=len(ws);ops=ws[i+1:i+count]
  if opcode==5:names[ops[0]]=struct.pack('<'+'I'*(len(ops)-1),*ops[1:]).split(b'\0')[0].decode()
  if opcode==15:entry=ops[1]
  if opcode==54:current=ops[1];functions[current]={'own':0,'calls':collections.Counter()}
  if current is not None:
   functions[current]['own']+=1
   if opcode==57:functions[current]['calls'][ops[2]]+=1;calls[(current,ops[2])]+=1
  if opcode==56:current=None
  i+=count
 @functools.cache
 def expanded(f):return functions[f]['own']+sum(n*expanded(c) for c,n in functions[f]['calls'].items())
 return {'function_count':len(functions),'call_sites':sum(calls.values()),'syntactic_expanded_entry_instructions':expanded(entry),
  'count_scope':'Function/body scaffolding and call instructions included; runtime loops are counted once, not multiplied by trip count; not an actual NIR/LLVM or memory measurement.',
  'entry_name':names.get(entry),'entry_own_instructions':functions[entry]['own']}
stages=[];native_proofs=[]
for kind,build in [('Camera',module.build_camera_program),('Transport',module.build_transport_program),('Endpoint',module.build_endpoint_program),('Dense',module.build_dense_program),('Initialize',module.build_initialize_program),('RayCamera',module.build_ray_camera_program)]:
 p=build(parallel=True);prefix=[p['instructions'],p['registers']]
 if kind not in ('Camera','RayCamera'):prefix.append(len(p['outputs']))
 expected=prefix+p['outputs']+p['operations']+p['layer_offsets']
 for text in headers:assert words(text,'k'+kind+'Program')==expected
 rows=[int(re.search('k'+kind+r'RowWords = (\d+);',text)[1]) for text in headers];assert rows[0]==rows[1]
 terms=4 if kind in ('Camera','RayCamera') else 5
 stem='ray_camera' if kind=='RayCamera' else kind.lower()
 for suffix,label in [('', ''),('_fp64','Fp64')]:
  path=HERE/'modules'/f'retained_{stem}{suffix}.spv';preset=PRESET/path.name
  assert path.read_bytes()==preset.read_bytes()
  assert words(headers[0],'k'+kind+label+'Shader')==words(headers[1],'k'+kind+label+'Shader')
  native_proofs.append({'path':str(path.relative_to(HERE)),'preset_path':str(preset.relative_to(ROOT)),'bytes':path.stat().st_size,'sha256':sha(path),'byte_identical':True})
 modes=[]
 for suffix,wide,label in [('_portable',False,'Portable'),('_portable_fp64',True,'PortableFp64')]:
  path=HERE/'modules'/f'retained_{stem}{suffix}.spv'
  data=path.read_bytes();assert list(struct.unpack('<'+str(len(data)//4)+'I',data))==words(headers[0],'k'+kind+label+'Shader')
  assembly=subprocess.run(['spirv-dis',str(path)],capture_output=True,text=True,check=True).stdout
  assert re.findall(r'OpCapability (\S+)',assembly)==['Shader']
  assert set(re.findall(r'OpTypeInt (\d+) [01]',assembly))=={'32'}
  assert 'OpTypeFloat' not in assembly
  entry=re.search(r'OpEntryPoint GLCompute (%\S+)',assembly)[1]
  assert 'OpExecutionMode '+entry+' LocalSize 1 1 1' in assembly
  assert 'OpControlBarrier' not in assembly
  assert not re.search(r'OpExecutionMode\S* .* (?:Denorm|RoundingMode|SignedZeroInfNan)',assembly)
  subprocess.run(['spirv-val','--target-env','vulkan1.2',str(path)],check=True,capture_output=True)
  modes.append({'label':'portable_fp64_request' if wide else 'portable_fp32_request','wide_product_definition':wide,'path':str(path.relative_to(HERE)),
   'sha256':sha(path),'bytes':len(data),'capabilities':['Shader'],'integer_widths':[32],'float_widths':[],
   'execution_modes':['LocalSize1x1x1'],'validation':'spirv-val --target-env vulkan1.2 PASS','static_calls':call_counts(path)})
 assert modes[0]['sha256']==modes[1]['sha256']
 stages.append({'kind':kind,'instructions':p['instructions'],'registers':p['registers'],'layers':len(p['layer_offsets'])-1,'terms':terms,
  'prefix_instructions':p.get('prefix_instructions',0),'row_words':rows[0],'shared_memory_bytes':p['registers']*terms*4+4,
  'program_sha256':hashlib.sha256(struct.pack('<'+'I'*len(expected),*expected)).hexdigest(),'program_matches_preset':True,
  'public_mode_modules_identical':True,'modes':modes})
report={'scope':'Isolated exact three-step portable Sum32 body sharing, based on the serial ordered-loop candidate. One invocation owns each ray; unchanged workgroup-per-ray dispatch. Offline validation only, no numerical qualification.',
 'source_revision':baseline['source_revision'],'baseline_sha256':sha(HERE/'baseline.json'),'dispatches':0,'pass':True,
 'changed_copied_sources':changed,'candidate_retained_triple_sha256':sha(HERE/'src/sirius/kernels/retained_triple.slang'),'native_source_branches_exactly_restored':True,
 'original_sources_unchanged':True,'copied_sources_sha256':{n:sha(HERE/n) for n in baseline['original_source_sha256']},
 'native_module_proofs':native_proofs,'native_modules_byte_identical_count':len(native_proofs),
 'program_words_unchanged':True,'scratch_source_unchanged':True,'row_words_unchanged':True,
 'abi_scope':'Arithmetic/scalar/scratch/program owners match the reviewed serial candidate; portable pair Sum32 shares the exact three original scalar-add operations and negated operands. Exact program arrays, roots, row-word constants and per-ray group dispatch match production. Portable execution is one invocation with serial layer evaluation and no control barriers; all native modules remain byte-identical.',
 'stages':stages,'canonical_transport_static_calls':call_counts(PRESET/'retained_transport_portable.spv')}
(HERE/'stage-artifact-report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'pass':True,'native_modules_identical':len(native_proofs),'portable_pairs_identical':len(stages),'report_sha256':sha(HERE/'stage-artifact-report.json'),
 'portable_modules':[{ 'stage':s['kind'],'bytes':s['modes'][0]['bytes'],'static_calls':s['modes'][0]['static_calls']} for s in stages],
 'canonical_transport_static_calls':report['canonical_transport_static_calls']}))
