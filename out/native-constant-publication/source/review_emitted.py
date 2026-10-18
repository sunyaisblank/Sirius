"""Self-contained finite saved-file proof; no project imports or tool/provider execution.
Typed normalization is adapted from the archived validation-only emitted reviewer.
Only function-local result IDs, recursively typed static IDs, and debug/function
names/occurrences identify counterparts; instruction order/operands remain exact.
"""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import difflib, hashlib, json, re

R=Path('/home/astra/.project/Sirius')
W=R/'out/native-constant-compiler-source'
P=R/'out/native-constant-publication'
ARCH=R/'attestations/qualification/5b2df58b/packed-retained-scratch-rejected/files/out/packed-retained-scratch-screen'
BASE=ARCH/'posttools/disassemble'
PRIOR=R/'attestations/qualification/5b2df58b/native-transport-validation-rejected/files/out/native-transport-validation-publication/source/review_emitted.py'
TOK=re.compile(r'"(?:\\.|[^"\\])*"|[^\s]+')
REV='5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e'
NAMES=['retained_transport','retained_endpoint','retained_dense','retained_initialize','retained_dopri_phase']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def identity(p):
 b=Path(p).read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def pin(p):return {'path':str(p),**identity(p)}
def signature(v):return hashlib.sha256(json.dumps(v,separators=(',',':')).encode()).hexdigest()
def match_pin(p,row):
 assert identity(p)=={k:row[k] for k in ('bytes','sha256')},str(p)
 if 'resolved_path' in row:assert str(Path(p).resolve())==row['resolved_path'],str(p)

class Module:
    def __init__(self,path):
        self.path=path; self.text=path.read_text();self.alias={}
        self.defs=dict(re.findall(r'^\s*(%\S+) = (.+)$',self.text,re.M))
        self.names=dict(re.findall(r'OpName (%\S+) "([^"]*)"',self.text))
        self.bodies=dict(re.findall(r'^\s*(%\S+) = (OpFunction\b.*?OpFunctionEnd)',self.text,re.M|re.S))
        self.keys={};self.by_key={};self.rows={};self.kept={};self.local={};counts=Counter()
        for ident,body in self.bodies.items():
            name=self.names[ident];key=(name,counts[name]);counts[name]+=1
            self.keys[ident]=key;self.by_key[key]=ident
            rows=[ident+' = '+body.splitlines()[0]]+body.splitlines()[1:]
            self.rows[key]=[TOK.findall(s) for s in rows if s.strip()]
            selected=list(enumerate(self.rows[key]))
            self.kept[key]=selected
            for i,row in selected:
                if len(row)>1 and row[1]=='=':self.local[row[0]]=(key,len([v for v in self.local.values() if v[0]==key]))
        self.globals={i:e for i,e in self.defs.items() if e.startswith('OpVariable ') and e.split()[2]!='Function'}
        for i in self.globals:
            if i not in self.names:
                builtin=re.findall(r'OpDecorate '+re.escape(i)+r' BuiltIn (\S+)',self.text)
                assert len(builtin)==1 and builtin[0] in ('WorkgroupId','LocalInvocationIndex')
                self.names[i]='BuiltIn:'+builtin[0]
        assert len({self.names[i] for i in self.globals})==len(self.globals)
        self.static={i:e for i,e in self.defs.items() if i not in self.keys and i not in self.local and e.startswith(('OpType','OpConstant','OpSpecConstant','OpExtInstImport ','OpUndef '))}
        self.norm=lru_cache(None)(self._norm)

    def _norm(self,i):
        if i in self.alias:return self.norm(self.alias[i])
        if i in self.keys:return ('function',self.keys[i])
        if i in self.globals:return ('global',self.names[i])
        if i in self.local:return ('local',self.local[i])
        assert i in self.static,(i,self.defs.get(i))
        return tuple(self.norm(v) if v.startswith('%') else v for v in TOK.findall(self.static[i]))

    def normalized_rows(self,key):
        return [tuple(self.norm(t) if t.startswith('%') else t for t in row) for _,row in self.kept[key]]

    def globals_normal(self):
        return {self.names[i]:tuple(self.norm(t) if t.startswith('%') else t for t in TOK.findall(e)) for i,e in self.globals.items()}

    def semantic_header(self):
        values=[]
        for line in self.text.splitlines():
            row=TOK.findall(line.strip())
            if not row:continue
            if row[0]=='OpEntryPoint':
                values.append(tuple(row[:2])+tuple([self.norm(row[2]),row[3],tuple(sorted(self.norm(t) for t in row[4:]))]))
            elif row[0] in ['OpCapability','OpExtension','OpMemoryModel','OpExecutionMode','OpExecutionModeId','OpDecorate','OpMemberDecorate','OpDecorateId']:
                values.append(tuple(self.norm(t) if t.startswith('%') else t for t in row))
        return Counter(values)


# Complete concrete Transport constructor witness, not a permissive pattern.
CONSTRUCTOR='''%RTConstant = OpFunction %RetainedTriple None %500
%x_2 = OpFunctionParameter %float
%1599 = OpLabel
%r_2 = OpVariable %_ptr_Function_RetainedTriple Function
%1600 = OpVariable %_ptr_Function_float Function
%1601 = OpFunctionCall %bool %RPBounded %x_2 %uint_2071986176
%1602 = OpLogicalNot %bool %1601
OpSelectionMerge %1603 None
OpBranchConditional %1602 %1604 %1603
%1604 = OpLabel
%1605 = OpFunctionCall %RetainedTriple %RTInvalid
OpReturnValue %1605
%1603 = OpLabel
%1606 = OpAccessChain %_ptr_Function_float %r_2 %int_0
%1607 = OpBitcast %uint %x_2
%1608 = OpBitwiseAnd %uint %1607 %uint_2147483647
%1609 = OpIEqual %bool %1608 %uint_0
OpSelectionMerge %1610 None
OpBranchConditional %1609 %1611 %1612
%1611 = OpLabel
%1613 = OpBitcast %float %uint_0
OpStore %1600 %1613
OpBranch %1610
%1612 = OpLabel
OpStore %1600 %x_2
OpBranch %1610
%1610 = OpLabel
%1614 = OpLoad %float %1600
OpStore %1606 %1614
%1615 = OpAccessChain %_ptr_Function_float %r_2 %int_1
%1616 = OpBitcast %float %uint_0
OpStore %1615 %1616
%1617 = OpAccessChain %_ptr_Function_float %r_2 %int_2
OpStore %1617 %1616
%1618 = OpAccessChain %_ptr_Function_float %r_2 %int_3
OpStore %1618 %1616
%1619 = OpAccessChain %_ptr_Function_uint %r_2 %int_4
OpStore %1619 %uint_1
%1620 = OpLoad %RetainedTriple %r_2
OpReturnValue %1620
OpFunctionEnd'''


# Fixed actual differences, not a pattern that silently removes unmatched code.
CSE={
 'retained_transport':{('ComputeMain',0):[(589,481)]},
 'retained_endpoint':{('Project',0):[(255,95),(256,95),(509,95)]},
 'retained_dense':{('FlatDense',0):[(158,134),(162,126),(336,272)]},
 'retained_initialize':{},'retained_dopri_phase':{}}

def cfg(m,key):
 rows=m.rows[key];blocks={};at=None;locations={}
 for i,row in enumerate(rows):
  op=row[2] if len(row)>1 and row[1]=='=' else row[0]
  if op=='OpLabel':at=row[0];blocks[at]=[]
  if at is not None and op!='OpFunctionEnd':blocks[at].append((i,row));locations[i]=at
 assert blocks
 successors={}
 for label,body in blocks.items():
  last=body[-1][1]
  if last[0]=='OpBranch':targets=last[1:]
  elif last[0]=='OpBranchConditional':targets=last[2:4]
  elif last[0]=='OpSwitch':targets=[last[2]]+last[4::2]
  else:assert last[0] in ('OpReturn','OpReturnValue','OpUnreachable');targets=[]
  assert set(targets)<=set(blocks);successors[label]=set(targets)
 entry=next(iter(blocks));reached={entry};todo=[entry]
 while todo:
  for label in successors[todo.pop()]-reached:reached.add(label);todo.append(label)
 assert reached==set(blocks)
 pred={label:set() for label in blocks}
 for label,targets in successors.items():
  for target in targets:pred[target].add(label)
 dom={label:({entry} if label==entry else set(blocks)) for label in blocks}
 changed=True
 while changed:
  changed=False
  for label in blocks:
   if label==entry:continue
   assert pred[label]
   value={label}|set.intersection(*(dom[x] for x in pred[label]))
   if value!=dom[label]:dom[label]=value;changed=True
 return locations,dom,pred

def apply_cse(a,b,name):
 witnesses=[]
 for key,pairs in CSE[name].items():
  rows=a.rows[key];locations,dom,pred=cfg(a,key);drops={d for d,k in pairs}
  assert len(drops)==len(pairs)
  for drop,keep in pairs:
   d=rows[drop];k=rows[keep]
   assert drop>keep and d[2:]==k[2:] and d[2]=='OpFunctionCall' and a.keys[d[4]]==('RTConstant',0)
   assert locations[keep] in dom[locations[drop]]
   if locations[keep]==locations[drop]:assert keep<drop
   a.alias[d[0]]=k[0]
   witnesses.append({'function':list(key),'removed_row':drop,'retained_row':keep,'removed_result':d[0],'retained_result':k[0],'identical_type_callee_argument_IDs':d[3:],'argument_definitions':[a.defs[v] for v in d[5:]],'removed_block':locations[drop],'retained_dominating_block':locations[keep],'removed_block_predecessors':sorted(pred[locations[drop]]),'retained_block_dominates_removed':True,'same_argument_SSA_values':True})
  a.kept[key]=[(i,row) for i,row in enumerate(rows) if i not in drops]
  for i in [i for i,v in a.local.items() if v[0]==key]:del a.local[i]
  count=0
  for i,row in a.kept[key]:
   if len(row)>1 and row[1]=='=':a.local[row[0]]=(key,count);count+=1
 a.norm.cache_clear()
 # The residual exact comparison below checks every substituted use as well.
 return witnesses

def pure_constructor(m):
 todo=[m.by_key[('RTConstant',0)]];visited=set();stores=loads=0
 allowed={'OpFunction','OpFunctionParameter','OpLabel','OpVariable','OpFunctionCall','OpLogicalNot','OpSelectionMerge','OpBranchConditional','OpReturnValue','OpAccessChain','OpBitcast','OpBitwiseAnd','OpIEqual','OpStore','OpBranch','OpLoad','OpFunctionEnd','OpULessThanEqual'}
 def function_pointer(ident):
  e=TOK.findall(m.defs[ident]);assert e[0] in ('OpVariable','OpAccessChain')
  pointer=TOK.findall(m.defs[e[1]]);assert pointer[:2]==['OpTypePointer','Function']
 while todo:
  ident=todo.pop()
  if ident in visited:continue
  visited.add(ident)
  for row in m.rows[m.keys[ident]]:
   op=row[2] if len(row)>1 and row[1]=='=' else row[0]
   assert op in allowed,('constructor_effect',op)
   if op=='OpFunctionParameter':assert m.defs[row[3]] in ('OpTypeFloat 32','OpTypeInt 32 0')
   if op=='OpFunctionCall':
    assert m.keys[row[4]] in [('RPBounded',0),('RTInvalid',0)];todo.append(row[4])
   if op=='OpStore':function_pointer(row[1]);stores+=1
   if op=='OpLoad':function_pointer(row[4]);loads+=1
   if op=='OpVariable':assert row[4]=='Function'
 assert {m.keys[i] for i in visited}=={('RTConstant',0),('RPBounded',0),('RTInvalid',0)}
 return {'closed_value_parameter_calltree':sorted([list(m.keys[i]) for i in visited]),'Function_only_stores':stores,'Function_only_loads':loads,'global_reads_or_writes':0,'atomics_or_barriers_or_extinst':0,'deterministic_bit_integer_and_control_operations_only':True,'initialized_loads':'Full exact constructor template writes temporary in both zero/nonzero arms and all5struct fields before load; unchanged RTInvalid writes all5fields before load.'}


def main():
 target=P/'independent-emitted-and-compiler-review.json';assert not target.exists()
 q=W/'candidate';pf=read(W/'root-preflight.json');before=read(q/'inputs-before.json');after=read(q/'inputs-after.json')
 assert len(pf['files'])==pf['static_input_count']==4450
 assert len(before)==pf['prospective_owner_input_count']==4451 and before==after
 for key,row in before.items():match_pin(Path(key) if Path(key).is_absolute() else R/key,row)
 for key,row in pf['files'].items():match_pin(Path(key),row)
 assert {str((Path(k) if Path(k).is_absolute() else R/k).resolve()) for k in before}=={str(Path(k).resolve()) for k in pf['files']}|{str((W/'root-preflight.json').resolve())}
 owner=read(q/'owner.json');assert owner['passed'] is True and owner['status']=='terminal' and owner['child_exit']==0 and owner['stop_reason'] is None
 assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[] and owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True
 assert owner['source']=={'revision':REV,'status':''} and owner['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
 assert owner['provider_mappings']==[] and owner['declared_output_environment']=={'GTEST_OUTPUT':''}
 assert owner['peak_sampled_rss_kib']<4194304 and owner['minimum_sampled_available_kib']>2097152 and owner['wall_seconds']<300
 assert len(owner['observed_births'])==7 and all(b['absent'] is True and b['current_start_ticks'] is None for b in owner['observed_births'])
 action=read(W/'actions.json')['candidate'];assert action['kind']=='build' and action['software_provider'] is False and action['expected_exit']==0
 assert action['limits']==owner['limits'] and owner['argv']==[a.replace('{folder}',str(q)) for a in action['argv']]
 recipe=read(W/'recipes.json');products=read(q/'products.json');results=read(q/'command-results.json')
 assert [p['name'] for p in recipe['profiles']]==NAMES
 assert products['schema']=='sirius-native-constant-five-compiler-products-v1' and products['passed'] is True
 assert products['runtime_admitted'] is False and products['emitted_semantics_admitted'] is False and products['optimizer_invoked'] is False
 assert len(results)==products['actual_tool_calls']==25 and set(products['products'])==set(NAMES)
 assert products['recipe']==pin(W/'recipes.json')
 oldproducts=read(ARCH/'root-compiler-products.json');assert oldproducts['passed'] is True
 olddis=read(BASE/'results.json');oldbefore=read(BASE/'inputs-before.json');oldafter=read(BASE/'inputs-after.json')
 oldtools=read(ARCH/'posttools/root-preflight.json')['files']
 mapping={};reports=[];diffs=[];template=None
 for index,profile in enumerate(recipe['profiles']):
  name=profile['name'];paths={s:q/(name+s) for s in recipe['expected_products_per_profile']};tools=recipe['tools'];src=Path(recipe['source_directory'])
  expected=[('compile',[tools['compiler'],str(src/(name+'.slang')),*profile['definitions'],'-I',str(src),'-O0','-target','spirv','-profile','spirv_1_5','-entry','ComputeMain','-stage','compute','-denorm-mode-fp32','preserve','-o',str(paths['.compiler.spv'])]),('raw-disassemble',[tools['disassembler'],str(paths['.compiler.spv']),'-o',str(paths['.raw.spvasm'])]),('assemble',[tools['assembler'],'--target-env','spv1.5',str(paths['.spvasm']),'-o',str(paths['.spv'])]),('validate',[tools['validator'],'--target-env','vulkan1.2',str(paths['.spv'])]),('final-disassemble',[tools['disassembler'],str(paths['.spv']),'-o',str(paths['.final.spvasm'])])]
  for actual,(label,argv) in zip(results[index*5:index*5+5],expected):
   assert actual['label']==name+'-'+label and actual['argv']==argv and actual['exit_code']==0
   for stream in ('stdout','stderr'):match_pin(Path(actual[stream]['path']),actual[stream])
  for suffix,path in paths.items():assert products['products'][name][suffix]==pin(path)
  annotated=paths['.spvasm'].read_text();inverse=re.sub(r'^.*OpCapability RoundingModeRTE\n','',annotated,flags=re.M);inverse=re.sub(r'^.*OpExecutionMode %\S+ RoundingModeRTE 32\n','',inverse,flags=re.M)
  assert inverse==paths['.raw.spvasm'].read_text()
  # Same saved disassembler and exact source SPV were used for canonical ASM.
  oldcmd,=[r for r in olddis if r.get('label')=='baseline-'+name]
  assert oldcmd['returncode']==0 if 'returncode' in oldcmd else oldcmd['exit_code']==0
  assert oldcmd['argv'][0]==tools['disassembler']
  match_pin(Path(tools['disassembler']),oldtools[tools['disassembler']])
  oldspv=oldcmd['argv'][1];oldkey=oldspv if oldspv in oldbefore else str(Path(oldspv).relative_to(R));assert oldbefore[oldkey]==oldafter[oldkey]
  baseline=pf['baseline_SPVs'][name+'.spv'];match_pin(Path(baseline['path']),oldbefore[oldkey])
  oldrow,=[r for r in oldproducts['products'] if r['name']==name+'.spv'];assert oldrow['baseline']==identity(Path(baseline['path']))
  basepath=BASE/('baseline-'+name+'.spvasm');assert oldcmd['argv'][-1].endswith('/baseline-'+name+'.spvasm')
  # Old saved command-result output identity, if present, must bind the ASM too.
  for field in ('output','assembly'):
   if field in oldcmd:match_pin(basepath,oldcmd[field])
  a=Module(basepath);b=Module(paths['.final.spvasm'])
  assert set(a.by_key)==set(b.by_key) and len(a.by_key)==len(b.by_key)
  assert a.globals_normal()==b.globals_normal(),('globals',name)
  assert Counter(a.norm(i) for i in a.static)==Counter(b.norm(i) for i in b.static),('static',name)
  # Header equality is checked after explicitly witnessed CSE remapping below.
  assert re.findall(r'OpExecutionMode %\S+ LocalSize (.+)',b.text)==['64 1 1']
  assert re.findall(r'OpExecutionMode %\S+ DenormPreserve (.+)',b.text)==['32'] and re.findall(r'OpExecutionMode %\S+ RoundingModeRTE (.+)',b.text)==['32']
  assert 'SignedZeroInfNanPreserve' not in b.text
  cse=apply_cse(a,b,name)
  assert a.semantic_header()==b.semantic_header(),('header',name)
  purity=pure_constructor(b)
  functions=[];hashes={};localpairs=[]
  for key in sorted(a.by_key):
   if key==('RTConstant',0):continue
   ar=a.normalized_rows(key);br=b.normalized_rows(key)
   assert ar==br,('body',name,key)
   ad=[row[0] for _,row in a.kept[key] if len(row)>1 and row[1]=='='];bd=[row[0] for _,row in b.kept[key] if len(row)>1 and row[1]=='=']
   assert len(ad)==len(bd);localpairs.append({'function':list(key),'baseline_id':a.by_key[key],'candidate_id':b.by_key[key],'definition_id_pairs':list(zip(ad,bd))})
   functions.append({'function':list(key),'rows':len(ar),'normalized_sha256':signature(ar)});hashes[str(key)]=signature(ar)
  assert len([k for k in a.by_key if k[0]=='RTConstant'])==1
  oldbody=a.rows[('RTConstant',0)];assert len(oldbody)==7 and oldbody[4][2]=='OpFunctionCall' and a.keys[oldbody[4][4]]==('RTPack',0)
  assert oldbody[4][5]==oldbody[1][0] and oldbody[4][6:]==[oldbody[3][0]]*3
  assert oldbody[3][2:] == ['OpBitcast','%float','%uint_0'] and a.defs['%uint_0']=='OpConstant %uint 0'
  if template is None:
   assert b.rows[('RTConstant',0)]==[TOK.findall(line) for line in CONSTRUCTOR.splitlines()]
   assert b.defs['%float']=='OpTypeFloat 32' and b.defs['%uint']=='OpTypeInt 32 0' and b.defs['%int']=='OpTypeInt 32 1'
   assert b.defs['%RetainedTriple']=='OpTypeStruct %float %float %float %float %uint'
   for ident,value in [('uint_0',0),('uint_1',1),('uint_2071986176',0x7b800000),('uint_2147483647',0x7fffffff)]:assert b.defs['%'+ident]==f'OpConstant %uint {value}'
   for value in range(5):assert b.defs['%int_'+str(value)]==f'OpConstant %int {value}'
   template=b.normalized_rows(('RTConstant',0))
  assert b.normalized_rows(('RTConstant',0))==template,('constructor',name)
  # Exact immutable globals/type/decorations and full outside bodies constitute ABI/workgroup/packing proof.
  staticpairs=[{'baseline_id':i,'candidate_ids':[j for j in b.static if b.norm(j)==a.norm(i)],'recursive_typed_signature':signature(a.norm(i))} for i in a.static]
  assert all(row['candidate_ids'] for row in staticpairs)
  mapping[name]={'functions':localpairs,'static':staticpairs,'globals':[{'name':a.names[i],'baseline_id':i,'candidate_id':next(j for j in b.globals if b.names[j]==a.names[i]),'typed_signature':signature(a.globals_normal()[a.names[i]])} for i in a.globals]}
  diffs.append(''.join(difflib.unified_diff([' '.join(x)+'\n' for x in oldbody],[' '.join(x)+'\n' for x in b.rows[('RTConstant',0)]],fromfile=name+'/canonical-RTConstant',tofile=name+'/candidate-RTConstant')))
  reports.append({'name':name,'baseline_SPV':baseline,'baseline_ASM':pin(basepath),'candidate_SPV':pin(paths['.spv']),'candidate_ASM':pin(paths['.final.spvasm']),'function_instances':len(a.by_key),'outside_RTConstant_instances_exact':len(functions),'ComputeMain_baseline_rows':len(a.rows[('ComputeMain',0)]),'ComputeMain_candidate_rows':len(b.rows[('ComputeMain',0)]),'ComputeMain_strict_alpha_exact':not any(c['function']==['ComputeMain',0] for c in cse),'explicit_CSE':cse,'constructor_transitive_purity':purity,'outside_instruction_rows_exact':sum(f['rows'] for f in functions),'outside_function_proof_digest':signature(functions),'constructor_rows':len(template),'typed_constructor_proof_digest':signature(template),'all_recursive_static_definitions_and_global_declarations_equal':True,'whole_capability_mode_decorations_multiset_equal':True,'local_size':[64,1,1]})
 idmap=P/'source/emitted-id-map.json';idmap.write_text(json.dumps(mapping,separators=(',',':'))+'\n')
 diff=P/'source/RTConstant.diff';diff.write_text('\n'.join(diffs))
 result={'schema':'sirius-native-constant-independent-emitted-and-compiler-review-v1','passed':True,'source_revision':REV,'reviewer':'/root/sanitizer_gate_contract','compiler_source_authored_by_reviewer':False,'prior_refinement_review_authored_by_reviewer':True,'actual_tools_operated_by_ROOT_only':True,
 'outside_RTConstant_functions_alpha_exact':False,'outside_constructor_semantics_preserved_with_explicit_pure_CSE':True,'public_ABI_equal':True,'bounded_native_constructor_refinement_validated':True,
 'compiler_owner':pin(q/'owner.json'),'compiler_products':pin(q/'products.json'),'command_results':pin(q/'command-results.json'),'modules':{n+'.spv':pin(q/(n+'.spv')) for n in NAMES},
 'terminal':{'static_current':4450,'before_after_current':4451,'actual_tool_calls':25,'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'minimum_sampled_available_kib':owner['minimum_sampled_available_kib'],'observed_births_recorded_absent':owner['observed_births'],'fresh_process_census_performed':False,'source_and_inputs_unchanged':True,'cleanup_errors':[],'remaining_owned':[]},
 'proof':reports,
 'constructor_contract':'Entire actual typed constructor matched across all5modules. First call RPBounded(x,0x7b800000); false returns unchanged RTInvalid. Surviving path uses absbits mask0x7fffffff and zero comparison; high selects bitcastuint0 for either zero sign, else exact x. Local struct members0/1/2/3/4 receive high/+0/+0/+0/uint1; complete struct loaded and returned. Only RPBounded and RTInvalid calls; old RTPack call removed. Bound/helper/invalid implementation unchanged. All alternatives/merges/stores/returns included; no extra constructor instructions ignored.',
 'normalization':'Saved strict-alpha proof/failure are preserved. Exactly7 explicitly enumerated duplicate RTConstant calls are removed (Transport1, Endpoint3, Dense3), each with identical raw argument IDs and a dominating retained call proven from the full closed CFG. Their result uses are substituted only by that retained SSA value. Every residual ordered instruction in every outside function, including ComputeMain, is then exact under typed alpha mapping. Function identity=(OpName, occurrence), local IDs paired by definition order. Types/constants/imports recurse through exact opcode/literal/operand signatures; global variables matched by names and separately exact typed declarations. Static signatures and semantic decoration/header multiplicities exactly equal. Only entry-interface declaration-list order treated as unordered same-global list; descriptor/layout operand order is preserved. OpName/debug order and assembler Bound are excluded. No arithmetic opcode, operand, branch, call, storage, barrier or decoration dropped.',
 'pins':{k:pin(p) for k,p in {'reader':Path(__file__),'prior_reader':PRIOR,'mapping':idmap,'readable_delta':diff,'before':q/'inputs-before.json','after':q/'inputs-after.json','preflight':W/'root-preflight.json','recipe':W/'recipes.json','source_refinement':R/'out/native-constant-source-review/semantic-refinement-review.json','strict_reader':P/'source/review_emitted.strict-alpha.py','strict_failure':P/'source/strict-alpha-finding.json','baseline_product_proof':ARCH/'root-compiler-products.json','baseline_disassembly_results':BASE/'results.json','baseline_disassembly_before':BASE/'inputs-before.json','baseline_disassembly_after':BASE/'inputs-after.json','baseline_disassembly_preflight':ARCH/'posttools/root-preflight.json'}.items()},
 'scope':'Five actual offline ordinary nativeFP32 modules only. Public modes/local64/storage/layout/packing/RMW/barriers/math/tables are unchanged. Full outside code is exact after only the7 explicitly proved pure-call CSE substitutions; strict alpha equality is false, not hidden. Constructor chooses a legal bounded canonical representation; universal prior-driver word equality is not claimed. Original same-provider32fullbuffer gate/all8controls/sole22benefit gates remain required. No provider, numerical, performance, adoption or qualification admission.',
 'runtime_admitted':False,'production_adopted':False,'operations':'Only passive standard-library saved-file/hash/JSON/text comparison; no project import/helper/compiler/SPIR-V tool/provider/test/network/Git mutation. Recorded owner birth-absence retained; no fresh process query.'}
 target.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'report':pin(target),'reader':pin(Path(__file__)),'mapping':pin(idmap),'readable_delta':pin(diff)},indent=2))

if __name__=='__main__':main()
