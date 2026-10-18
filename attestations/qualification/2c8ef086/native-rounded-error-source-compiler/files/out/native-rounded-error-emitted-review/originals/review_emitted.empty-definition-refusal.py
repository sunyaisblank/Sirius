"""Passive exact saved-SPIRV comparison; standard library only, no project imports.
Module normalization is preserved from the fd530 strict reader. No pure-CSE or
RTConstant exemption is imported. ROOT owns actual compiler and SPIR-V tools.
"""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import hashlib,json,re,difflib
R=Path(__file__).resolve().parents[2]
W=R/'out/native-rounded-error-compiler-source'
P=Path(__file__).resolve().parent
BASE=R/'attestations/qualification/5b2df58b/native-constant-controls-benefit-pass/files/out/native-constant-compiler-source/candidate'
TOK=re.compile(r'"(?:\\.|[^"\\])*"|[^\s]+')
NAMES=['retained_transport','retained_endpoint','retained_dense','retained_initialize','retained_dopri_phase']
def identity(p):
 b=Path(p).read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def pin(p):return {'path':str(p),**identity(p)}
def signature(value):return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))

class Module:
    def __init__(self,path):
        self.path=path; self.text=path.read_text()
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



GUARD_TEMPLATE="""$MAG = OpBitcast @uint $ABS
$LOW = OpUGreaterThanEqual @bool $MAG @u234881024
OpSelectionMerge $GJOIN None
OpBranchConditional $LOW $LOWOK $LOWFAIL
$LOWOK = OpLabel
$UPPER = OpULessThanEqual @bool $MAG @u1585446912
OpStore $BOOL $UPPER
OpBranch $GJOIN
$LOWFAIL = OpLabel
OpStore $BOOL @false
OpBranch $GJOIN
$GJOIN = OpLabel
$GUARD = OpLoad @bool $BOOL
OpSelectionMerge $EJOIN None
OpBranchConditional $GUARD $TAKE $FALL
$TAKE = OpLabel
$SUB = OpISub @uint $MAG @u192937984
$ADD = OpIAdd @uint $SUB @u2
$TRANS = OpBitcast @float $ADD
OpStore $ERROR $TRANS
OpBranch $EJOIN
$FALL = OpLabel
$SCALE = OpBitcast @float @u872415232
$MUL = OpFunctionCall @float @fnRPUpMultiply $ABS $SCALE
$TINY = OpBitcast @float @u5
$FALLERROR = OpFunctionCall @float @fnRPUpAdd $MUL $TINY
OpStore $ERROR $FALLERROR
OpBranch $EJOIN
$EJOIN = OpLabel"""


def expected_static(marker):
 if marker in ('@uint','@float','@bool'):
  return {'@uint':('OpTypeInt','32','0'),'@float':('OpTypeFloat','32'),'@bool':('OpTypeBool',)}[marker]
 if marker=='@false':return ('OpConstantFalse',('OpTypeBool',))
 if marker.startswith('@u'):return ('OpConstant',('OpTypeInt','32','0'),marker[2:])
 if marker.startswith('@fn'):return ('function',(marker[3:],0))
 raise ValueError(marker)


def match_guard(m,rows,roles):
 expected=[line.split() for line in GUARD_TEMPLATE.splitlines()]
 assert len(rows)==len(expected)
 for row,template in zip(rows,expected):
  assert len(row)==len(template),(row,template)
  for token,want in zip(row,template):
   if want.startswith('$'):
    if want not in roles:
     assert token.startswith('%') and token not in roles.values(),(want,token)
     roles[want]=token
    assert token==roles[want],(want,token,roles[want])
   elif want.startswith('@'):
    assert m.norm(token)==expected_static(want),(want,m.norm(token))
   else:assert token==want,(token,want)
 return roles


def function_pointer_type(m,ident,kind):
 row=TOK.findall(m.defs[ident]);assert row[0]=='OpVariable' and row[2:] == ['Function']
 assert m.norm(row[1])==('OpTypePointer','Function',expected_static(kind))


def validate_and_restore_helper(a,b,key):
 rows=b.rows[key]
 abs_index,= [i for i,row in enumerate(rows) if len(row)>4 and row[2]=='OpFunctionCall' and b.keys.get(row[4])==('RPAbs',0)]
 term_index,= [i for i,row in enumerate(rows) if len(row)>4 and row[2]=='OpFunctionCall' and b.keys.get(row[4],('',0))[0]=='RTAccumulateTerm']
 error_id,= [row[0] for row in rows if len(row)>4 and row[2]=='OpVariable' and
            b.norm(row[3])==('OpTypePointer','Function',('OpTypeFloat','32'))]
 bool_id,= [row[0] for row in rows if len(row)>4 and row[2]=='OpVariable' and
            b.norm(row[3])==('OpTypePointer','Function',('OpTypeBool',))]
 function_pointer_type(b,error_id,'@float');function_pointer_type(b,bool_id,'@bool')
 roles=match_guard(b,rows[abs_index+1:term_index],{'$ABS':rows[abs_index][0],'$ERROR':error_id,'$BOOL':bool_id})
 assert all(bool_id not in row for row in rows[term_index:]),'guard bool escaped into tail'
 error_loads=[(i,row) for i,row in enumerate(rows) if len(row)>4 and row[2]=='OpLoad' and row[4]==error_id]
 assert len(error_loads)==1 and error_loads[0][0]>term_index
 load_index,load_row=error_loads[0]
 occurrences=[(i,row) for i,row in enumerate(rows) if error_id in row]
 assert len(occurrences)==4 and occurrences[0][1][0]==error_id and sum(row[0]=='OpStore' for _,row in occurrences)==2
 # All newly emitted instructions match the exact typed guard template above.
 # Substitute only the source-proved equal error value, restoring the original
 # fallback's SSA expression for a complete instruction-order comparison.
 kept=[]
 retain={roles[x] for x in ('$SCALE','$MUL','$TINY','$FALLERROR')}
 for i,row in enumerate(rows):
  if row[0]==error_id or i==load_index:continue
  if abs_index<i<term_index and row[0] not in retain:continue
  kept.append((i,[roles['$FALLERROR'] if v==load_row[0] else v for v in row]))
 b.kept[key]=kept
 for ident,location in list(b.local.items()):
  if location[0]==key:del b.local[ident]
 counter=0
 for _,row in kept:
  if len(row)>1 and row[1]=='=':b.local[row[0]]=(key,counter);counter+=1
 b.norm.cache_clear()
 assert a.normalized_rows(key)==b.normalized_rows(key),('helper inverse',key)
 return {'function':list(key),'canonical_rows':len(a.rows[key]),'candidate_rows':len(rows),
         'exact_guard_rows':len(GUARD_TEMPLATE.splitlines()),'roles':roles,
         'restored_all_original_rows_and_uses_exact':True}


def match(p,row):assert identity(p)=={k:row[k] for k in ('bytes','sha256')},str(p)


def main():
 target=P/'independent-emitted-review.json';assert not target.exists()
 q=W/'candidate';pf=read(W/'root-preflight.json');owner=read(q/'owner.json')
 before=read(q/'inputs-before.json');after=read(q/'inputs-after.json')
 assert before==after and len(before)==pf['prospective_owner_input_count']
 assert owner['passed'] is True and owner['status']=='terminal' and owner['child_exit']==0 and owner['stop_reason'] is None
 assert owner['source']=={'revision':'2c8ef08693d46f57d5dea14efdbddd4874dac684','status':''}
 assert owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True
 assert owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[] and owner['provider_mappings']==[]
 assert owner['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152}
 assert owner['wall_seconds']<300 and owner['peak_sampled_rss_kib']<4194304 and owner['minimum_sampled_available_kib']>2097152
 assert all(row['absent'] is True and row['current_start_ticks'] is None for row in owner['observed_births'])
 products=read(q/'products.json');commands=read(q/'command-results.json');recipe=read(W/'recipes.json')
 assert products['schema']=='sirius-native-rounded-error-five-compiler-products-v1' and products['passed'] is True
 assert products['runtime_admitted'] is False and products['emitted_semantics_admitted'] is False and products['optimizer_invoked'] is False
 assert len(commands)==products['actual_tool_calls']==25 and set(products['products'])==set(NAMES)
 assert products['recipe']==pin(W/'recipes.json')
 oldproducts=read(BASE/'products.json');oldcommands=read(BASE/'command-results.json');oldbefore=read(BASE/'inputs-before.json')
 baseline_products=read(W/'baseline-products.json')
 source=W/'source/src/sirius/kernels';source_files=sorted(p for p in source.iterdir() if p.is_file())
 assert len(source_files)==33
 for path in source_files:
  expected=R/'out/native-rounded-error-source/retained_triple.slang' if path.name=='retained_triple.slang' else R/'src/sirius/kernels'/path.name
  assert path.read_bytes()==expected.read_bytes(),('private source',path.name)
 reports=[];helper_proofs=[];deltas=[]
 for index,name in enumerate(NAMES):
  profile=recipe['profiles'][index];assert profile['name']==name and profile['definitions']==[]
  paths={suffix:q/(name+suffix) for suffix in recipe['expected_products_per_profile']};tools=recipe['tools']
  expected=[('compile',[tools['compiler'],str(source/(name+'.slang')),*profile['definitions'],'-I',str(source),'-O0','-target','spirv','-profile','spirv_1_5','-entry','ComputeMain','-stage','compute','-denorm-mode-fp32','preserve','-o',str(paths['.compiler.spv'])]),('raw-disassemble',[tools['disassembler'],str(paths['.compiler.spv']),'-o',str(paths['.raw.spvasm'])]),('assemble',[tools['assembler'],'--target-env','spv1.5',str(paths['.spvasm']),'-o',str(paths['.spv'])]),('validate',[tools['validator'],'--target-env','vulkan1.2',str(paths['.spv'])]),('final-disassemble',[tools['disassembler'],str(paths['.spv']),'-o',str(paths['.final.spvasm'])])]
  for row,(label,argv) in zip(commands[index*5:index*5+5],expected):
   assert row['label']==name+'-'+label and row['argv']==argv and row['exit_code']==0
   for stream in ('stdout','stderr'):match(Path(row[stream]['path']),row[stream])
  for suffix,path in paths.items():assert products['products'][name][suffix]==pin(path)
  annotated=paths['.spvasm'].read_text()
  inverse=re.sub(r'^.*OpCapability RoundingModeRTE\n','',annotated,flags=re.M)
  inverse=re.sub(r'^.*OpExecutionMode %\S+ RoundingModeRTE 32\n','',inverse,flags=re.M)
  assert inverse==paths['.raw.spvasm'].read_text()
  baseline=BASE/(name+'.spv');baseasm=BASE/(name+'.final.spvasm')
  match(baseline,oldproducts['products'][name]['.spv']);match(baseasm,oldproducts['products'][name]['.final.spvasm'])
  canonical=R/baseline_products['SPVs'][name+'.spv']['path'];match(canonical,baseline_products['SPVs'][name+'.spv'])
  assert baseline.read_bytes()==canonical.read_bytes()
  oldcmd,=[row for row in oldcommands if row['label']==name+'-final-disassemble']
  assert oldcmd['exit_code']==0 and oldcmd['argv'][0]==tools['disassembler']
  oldtool=oldbefore[tools['disassembler']];match(Path(tools['disassembler']),oldtool)
  a,b=Module(baseasm),Module(paths['.final.spvasm'])
  assert set(a.by_key)==set(b.by_key) and a.globals_normal()==b.globals_normal()
  extra=Counter(b.norm(i) for i in b.static)-Counter(a.norm(i) for i in a.static)
  assert not Counter(a.norm(i) for i in a.static)-Counter(b.norm(i) for i in b.static)
  assert extra==Counter({('OpConstant',('OpTypeInt','32','0'),'234881024'):1,('OpConstant',('OpTypeInt','32','0'),'192937984'):1}),('extra static',name,extra)
  outside=[]
  for key in sorted(a.by_key):
   if key[0]=='RTAccumulateRoundedProduct':continue
   assert a.normalized_rows(key)==b.normalized_rows(key),('outside strict alpha',name,key)
   outside.append({'function':list(key),'rows':len(a.rows[key]),'digest':signature(a.normalized_rows(key))})
  helpers=[key for key in a.by_key if key[0]=='RTAccumulateRoundedProduct'];assert helpers
  for key in helpers:
   deltas.append(''.join(difflib.unified_diff([' '.join(row)+'\n' for row in a.rows[key]],[' '.join(row)+'\n' for row in b.rows[key]],fromfile=name+'/canonical-'+str(key),tofile=name+'/candidate-'+str(key))))
   helper_proofs.append({'module':name,**validate_and_restore_helper(a,b,key)})
  # After the exact helper inverse, result-local NoContraction IDs have their
  # original roles. No decorations are omitted from this complete comparison.
  assert a.semantic_header()==b.semantic_header(),('capability/ABI/decorations',name)
  for op in ('FAdd','FSub','FMul'):
   ids=re.findall(r'(%\S+) = Op'+op+r' ',b.text)
   decorations=set(re.findall(r'OpDecorate (%\S+) NoContraction',b.text))
   assert ids and set(ids)<=decorations
  assert re.findall(r'OpExecutionMode %\S+ LocalSize (.+)',b.text)==['64 1 1']
  assert re.findall(r'OpExecutionMode %\S+ DenormPreserve (.+)',b.text)==['32']
  assert re.findall(r'OpExecutionMode %\S+ RoundingModeRTE (.+)',b.text)==['32']
  assert not re.search(r'\b(?:FPFastMathMode|FPFastMathDefault|RelaxedPrecision|OpFmaKHR)\b',b.text)
  assert not re.search(r'OpType(?:Float|Int) 64',b.text)
  reports.append({'name':name,'canonical_SPV':pin(canonical),'saved_baseline_SPV':pin(baseline),'saved_baseline_ASM':pin(baseasm),'candidate_SPV':pin(paths['.spv']),'candidate_ASM':pin(paths['.final.spvasm']),'functions':len(a.by_key),'outside_helper_instances_exact':len(outside),'outside_instruction_rows_exact':sum(row['rows'] for row in outside),'outside_proof_digest':signature(outside),'helper_instances':len(helpers),'full_header_and_global_ABI_equal':True,'mode':'local64x1x1/RTE32/DenormPreserve32/NoContraction; no relaxed/wide/fused arithmetic'})
 diff=P/'guarded-helper.diff';diff.write_text('\n'.join(deltas))
 witnesses=P/'helper-typed-witnesses.json';witnesses.write_text(json.dumps(helper_proofs,indent=2)+'\n')
 result={'schema':'sirius-native-rounded-error-independent-emitted-review-v1','passed':True,'source_revision':'2c8ef08693d46f57d5dea14efdbddd4874dac684','outside_RTAccumulateRoundedProduct_functions_alpha_exact':True,'outside_RTAccumulateRoundedProduct_semantics_preserved':True,'public_ABI_equal':True,'bounded_native_rounded_product_error_equivalence_validated':True,'guarded_native_error_word_refinement_validated':True,'compiler_owner':pin(q/'owner.json'),'compiler_products':pin(q/'products.json'),'command_results':pin(q/'command-results.json'),'modules':{name+'.spv':pin(q/(name+'.spv')) for name in NAMES},'proof':reports,'helper_typed_witnesses':pin(witnesses),'exact_helper_delta':pin(diff),'reader':pin(Path(__file__)),'prior_strict_reader':pin(P/'originals/review_emitted.strict-alpha.py'),'source_equivalence_review':pin(R/'out/native-rounded-error-source/independent-source-review.json'),'root_compiler_terminal':pin(R/'out/native-rounded-error-publication/root-compiler-terminal.json'),'terminal_scope':{'before_after_rows_equal':len(before),'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'recorded_births_absent':len(owner['observed_births']),'new_process_query':False,'all_current_host_inputs_rehashed_by_this_reader':False},'normalization':'Original strict typed normalization reused. Every outside function including RTConstant and ComputeMain is strictly alpha-exact with no CSE/side-effect exception. Only two new uint constants occur. Every changed helper matches the entire exact29-row typed guard/word/fallback CFG template, with private float error local initialized by both arms before sole tail load. Removing only that proved insertion, replacing that one error load by original fallback SSA, and retaining its four exact original arithmetic instructions reproduces every original helper instruction and use in order. Complete capability/mode/type/global/layout/decorations multiset then equals canonical; no NoContraction or other decoration is dropped.','baseline_provenance':'Each saved fd530 candidate SPV/full final ASM is checked against its original products record and successful disassembler record; the saved SPV equals the actual current2c8 canonical module byte-for-byte. Selected disassembler identity matches original sealed tool. No old live raw path is demanded.','scope':'Five ordinary native modules only; exact source guarded error expression and typed emitted fallback/ABI/modes. Same-provider full buffers/original8 controls/sole original22 benefit remain required. No separate lowering allowance, inherited RTConstant refinement or CSE exception.','runtime_admitted':False,'production_adopted':False,'operations':'Passive standalone standard-library JSON/hash/ASM comparison only; no project import/helper/compiler/SPIR-V tool/provider/test/network/Git mutation. Existing ROOT terminal seals reused; no new all-host sweep.'}
 target.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'report':pin(target),'witnesses':pin(witnesses),'diff':pin(diff),'reader':pin(Path(__file__))}))

if __name__=='__main__':main()
