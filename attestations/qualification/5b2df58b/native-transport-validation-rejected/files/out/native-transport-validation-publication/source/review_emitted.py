"""Finite saved-text/data review only: no project imports, subprocess or GPU/tool call."""
from pathlib import Path
from collections import Counter, defaultdict
from functools import lru_cache
import difflib, hashlib, json, re

R=Path('/home/astra/.project/Sirius')
W=R/'out/native-transport-validation-compiler-source'
P=R/'out/native-transport-validation-publication'
BASE=R/'attestations/qualification/5b2df58b/packed-retained-scratch-rejected/files/out/packed-retained-scratch-screen/posttools/disassemble/baseline-retained_transport.spvasm'
TOK=re.compile(r'"(?:\\.|[^"\\])*"|[^\s]+')


def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def identity(p):
    b=Path(p).read_bytes()
    return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def pin(p): return {'path':str(p),**identity(p)}
def signature(v): return hashlib.sha256(json.dumps(v,separators=(',',':')).encode()).hexdigest()


class Module:
    def __init__(self,path,cut):
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
            selected=[(i,row) for i,row in enumerate(self.rows[key]) if key[0]!='ComputeMain' or not any(a<=i<b for a,b in cut)]
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


def check_validation(m):
    rows=m.rows[('ComputeMain',0)]
    reset='''%754 = OpLoad %uint %gl_LocalInvocationIndex
%755 = OpIEqual %bool %754 %uint_0
OpSelectionMerge %756 None
OpBranchConditional %755 %757 %756
%757 = OpLabel
OpStore %retainedStatus %uint_0
OpBranch %756
%756 = OpLabel'''
    changed='''%797 = OpLoad %uint %gl_LocalInvocationIndex
OpStore %i %797
OpBranch %798
%798 = OpLabel
OpLoopMerge %799 %800 None
OpBranch %801
%801 = OpLabel
OpBranch %802
%802 = OpLabel
OpSelectionMerge %803 None
OpSwitch %int_0 %804
%804 = OpLabel
%805 = OpLoad %uint %i
%806 = OpULessThan %bool %805 %uint_46
OpSelectionMerge %807 None
OpBranchConditional %806 %807 %808
%807 = OpLabel
%809 = OpLoad %uint %i
%810 = OpFunctionCall %RetainedTriple %ReadOriginal %809 %row
%811 = OpCompositeExtract %uint %810 4
%812 = OpIEqual %bool %811 %uint_0
OpSelectionMerge %813 None
OpBranchConditional %812 %814 %813
%814 = OpLabel
%815 = OpFunctionCall %void %RejectRetainedProgram
OpBranch %813
%813 = OpLabel
OpBranch %803
%803 = OpLabel
OpBranch %816
%816 = OpLabel
%817 = OpLoad %uint %i
%i_1 = OpIAdd %uint %817 %uint_64
OpStore %i %i_1
OpBranch %800
%800 = OpLabel
OpBranch %798
%808 = OpLabel
OpBranch %799
%799 = OpLabel
OpControlBarrier %uint_2 %uint_2 %uint_264
%818 = OpLoad %uint %retainedStatus
%819 = OpINotEqual %bool %818 %uint_0
OpSelectionMerge %820 None
OpBranchConditional %819 %821 %820
%821 = OpLabel
OpReturn
%820 = OpLabel'''
    assert rows[204:212]==[s.split() for s in reset.splitlines()]
    assert rows[285:333]==[s.split() for s in changed.splitlines()]
    assert rows[212]==['OpControlBarrier','%uint_2','%uint_1','%uint_2376']
    assert rows[333]==['%chart','=','OpFunctionCall','%RetainedTriple','%ReadOriginal','%uint_44','%row']
    for ident,value in [('%uint_0',0),('%uint_1',1),('%uint_2',2),('%uint_46',46),('%uint_64',64),('%uint_264',264),('%uint_2376',2376)]:assert m.defs[ident]==f'OpConstant %uint {value}'
    assert m.defs['%uint']=='OpTypeInt 32 0'
    assert m.defs['%retainedStatus']=='OpVariable %_ptr_Workgroup_uint Workgroup' and m.defs['%_ptr_Workgroup_uint']=='OpTypePointer Workgroup %uint'
    assert m.defs['%int_2295']=='OpConstant %int 2295' and m.defs['%_arr_uint_int_2295']=='OpTypeArray %uint %int_2295'
    assert m.defs['%retainedScratch']=='OpVariable %_ptr_Workgroup__arr_uint_int_2295 Workgroup'
    atomic=m.rows[('RejectRetainedProgram',0)]
    assert len(atomic)==5 and atomic[2][2:]==['OpAtomicOr','%uint','%retainedStatus','%uint_1','%uint_0','%uint_1'] and atomic[3:]==[['OpReturn'],['OpFunctionEnd']]
    todo=[m.by_key[('ReadOriginal',0)]];visited=set();stores=0
    while todo:
        ident=todo.pop()
        if ident in visited:continue
        visited.add(ident)
        for row in m.rows[m.keys[ident]]:
            op=row[2] if len(row)>1 and row[1]=='=' else row[0]
            assert not op.startswith(('OpAtomic','OpImageWrite','OpCopyMemory','OpControlBarrier','OpMemoryBarrier'))
            if op=='OpFunctionCall':todo.append(row[4])
            if op=='OpStore':
                target=row[1];expression=TOK.findall(m.defs[target]);pointer=expression[1]
                assert TOK.findall(m.defs[pointer])[:2]==['OpTypePointer','Function'],(target,m.defs[target])
                stores+=1
    indices=[lane for lane in range(64) if lane<46]
    assert indices==list(range(46)) and max(i+64 for i in indices)==109
    return reset,changed,{'ReadOriginal_transitive_function_instances':sorted([list(m.keys[i]) for i in visited]),'Function_storage_only_stores':stores,'coverage_indices':indices,'max_after_increment':109}


def main():
    target=P/'independent-emitted-and-compiler-review.json';assert not target.exists()
    before=read(W/'candidate/inputs-before.json');after=read(W/'candidate/inputs-after.json');pf=read(W/'root-preflight.json');assert len(before)==4489 and before==after and len(pf['files'])==4488
    for key,row in before.items():
        p=Path(key) if Path(key).is_absolute() else R/key
        assert identity(p)=={k:row[k] for k in ['bytes','sha256']} and str(p.resolve())==row['resolved_path'],key
    for key,row in pf['files'].items():assert identity(Path(key))=={k:row[k] for k in ['bytes','sha256']},key
    owner=read(W/'candidate/owner.json');assert owner['passed'] is True and owner['status']=='terminal' and owner['child_exit']==0 and owner['stop_reason'] is None and owner['cleanup_errors']==[] and owner['remaining_owned_processes']==[]
    assert owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True and owner['source']=={'revision':'5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e','status':''}
    assert owner['limits']=={'seconds':300,'rss_kib':4194304,'reserve_kib':2097152} and owner['provider_mappings']==[] and owner['declared_output_environment']=={'GTEST_OUTPUT':''}
    births=[]
    for b in owner['observed_births']:
        assert b['absent'] is True
        try:
            s=Path(f"/proc/{b['pid']}/stat").read_text();current=s[s.rfind(')')+2:].split()[19]
        except FileNotFoundError:current=None
        assert str(current)!=str(b['saved_start_ticks']);births.append({'pid':b['pid'],'saved_start_ticks':str(b['saved_start_ticks']),'fresh_current_start_ticks':current,'recorded_birth_absent':True})
    results=read(W/'candidate/command-results.json');products=read(W/'candidate/products.json');recipe=read(W/'recipes.json');assert len(results)==products['actual_tool_calls']==5 and products['passed'] is True and products['optimizer_invoked'] is False
    assert products['runtime_admitted'] is False and products['emitted_validation_semantics_admitted'] is False and products['schema']=='sirius-native-transport-validation-one-compiler-products-v1'
    profile,=recipe['profiles'];assert profile['name']=='retained_transport'
    q=W/'candidate';paths={s:q/('retained_transport'+s) for s in recipe['expected_products_per_profile']};tools=recipe['tools'];src=Path(recipe['source_directory'])
    expected=[('compile',[tools['compiler'],str(src/'retained_transport.slang'),*profile['definitions'],'-I',str(src),'-O0','-target','spirv','-profile','spirv_1_5','-entry','ComputeMain','-stage','compute','-denorm-mode-fp32','preserve','-o',str(paths['.compiler.spv'])]),('raw-disassemble',[tools['disassembler'],str(paths['.compiler.spv']),'-o',str(paths['.raw.spvasm'])]),('assemble',[tools['assembler'],'--target-env','spv1.5',str(paths['.spvasm']),'-o',str(paths['.spv'])]),('validate',[tools['validator'],'--target-env','vulkan1.2',str(paths['.spv'])]),('final-disassemble',[tools['disassembler'],str(paths['.spv']),'-o',str(paths['.final.spvasm'])])]
    for actual,(label,argv) in zip(results,expected):
        assert actual['label']=='retained_transport-'+label and actual['argv']==argv and actual['exit_code']==0
        for stream in ['stdout','stderr']:assert identity(Path(actual[stream]['path']))=={k:actual[stream][k] for k in ['bytes','sha256']}
    for suffix,p in paths.items():assert products['products']['retained_transport'][suffix]==pin(p)
    annotated=paths['.spvasm'].read_text();inverse=re.sub(r'^.*OpCapability RoundingModeRTE\n','',annotated,flags=re.M);inverse=re.sub(r'^.*OpExecutionMode %\S+ RoundingModeRTE 32\n','',inverse,flags=re.M);assert inverse==paths['.raw.spvasm'].read_text()
    a=Module(BASE,[(277,315)]);b=Module(paths['.final.spvasm'],[(204,212),(285,333)])
    assert len(a.by_key)==len(b.by_key)==70 and set(a.by_key)==set(b.by_key)
    reset,changed,purity=check_validation(b)
    assert a.globals_normal()==b.globals_normal()
    assert Counter(a.norm(i) for i in a.static)==Counter(b.norm(i) for i in b.static)
    assert a.semantic_header()==b.semantic_header()
    assert re.findall(r'OpExecutionMode %\S+ LocalSize (.+)',b.text)==['64 1 1']
    mapping={'functions':[],'global_variables':[],'static_type_constant_imports':[]};bodyproof=[]
    for key in sorted(a.by_key):
        old=a.normalized_rows(key);new=b.normalized_rows(key);assert old==new,('function',key)
        olddefs=[row[0] for _,row in a.kept[key] if len(row)>1 and row[1]=='='];newdefs=[row[0] for _,row in b.kept[key] if len(row)>1 and row[1]=='='];assert len(olddefs)==len(newdefs)
        mapping['functions'].append({'function':list(key),'old_id':a.by_key[key],'candidate_id':b.by_key[key],'local_definition_id_pairs':list(zip(olddefs,newdefs))})
        bodyproof.append({'function':list(key),'kept_instruction_rows':len(old),'normalized_sha256':signature(old)})
    assert len(a.kept[('ComputeMain',0)])==len(b.kept[('ComputeMain',0)])==1141
    for old in a.globals:
        other,=[i for i in b.globals if b.names[i]==a.names[old]];mapping['global_variables'].append({'name':a.names[old],'old_id':old,'candidate_id':other,'typed_declaration_signature':signature(a.globals_normal()[a.names[old]])})
    for ident in a.static:
        alternatives=[i for i in b.static if b.norm(i)==a.norm(ident)]
        assert alternatives
        mapping['static_type_constant_imports'].append({'old_id':ident,'candidate_ids':alternatives,'exact_recursive_signature_sha256':signature(a.norm(ident))})
    amain=a.rows[('ComputeMain',0)];bmain=b.rows[('ComputeMain',0)]
    diff=P/'source/validation-regions.diff';diff.write_text(''.join(difflib.unified_diff([' '.join(x)+'\n' for x in amain[277:315]],[' '.join(x)+'\n' for x in bmain[285:333]],fromfile='canonical-main-rows277-314',tofile='candidate-main-rows285-332'))+'\nADDED BEFORE UNCHANGED INITIAL ALL-MEMORY JOIN:\n'+reset+'\n')
    idmap=P/'source/emitted-id-map.json';idmap.write_text(json.dumps(mapping,indent=2)+'\n')
    result={'schema':'sirius-native-transport-validation-independent-emitted-and-compiler-review-v1','passed':True,'source_revision':owner['source']['revision'],'reviewer':'/root/sanitizer_gate_contract','compiler_source_authored_by_reviewer':True,'actual_tools_operated_by_ROOT_only':True,'scope':'Saved actual nativeTransport compiler/typed-IR comparison only; no runtime/numerical/performance/adoption admission.',
      'pins':{k:pin(p) for k,p in {'compiler_owner':q/'owner.json','compiler_products':q/'products.json','command_results':q/'command-results.json','before':q/'inputs-before.json','after':q/'inputs-after.json','preflight':W/'root-preflight.json','recipe':W/'recipes.json','baseline_SPV':W/'originals/canonical-retained_transport.spv','baseline_ASM':BASE,'candidate_SPV':paths['.spv'],'candidate_ASM':paths['.final.spvasm'],'reader':Path(__file__),'mapping':idmap,'readable_delta':diff}.items()},
      'terminal':{'static_current':4488,'owner_before_after_current':4489,'actual_tool_calls':5,'wall_seconds':owner['wall_seconds'],'peak_sampled_rss_kib':owner['peak_sampled_rss_kib'],'fresh_recorded_births':births,'cleanup_errors':[],'remaining_owned':[]},
      'unchanged_numerical_helpers_alpha_exact':True,'public_ABI_equal':True,'validation_only_CFG_validated':True,
      'proof':{'complete_function_instances':70,'helper_instances_unchanged':69,'ComputeMain_kept_rows_exact':1141,'removed_baseline_main_range':[277,315],'candidate_added_reset_range':[204,212],'candidate_replacement_validation_range':[285,333],'whole_globals_recursive_type_constant_import_multiset_equal':True,'all_capabilities_modes_bindings_layout_and_NoContraction_decorations_multiset_equal':True,'entry_interface':'Same8globals; retainedStatus/retainedScratch declaration-list order changes, which does not change entry ABI. No ordered descriptor/binding/layout operand is dropped or sorted.','static_uint264':'Already existed baseline; only declaration order moved. AcquireRelease|WorkgroupMemory barrier operand value unchanged.','local_size':[64,1,1],'scratch_uint_words':2295,'scratch_bytes':9180,'body_signatures':bodyproof,'ReadOriginal_pure_call_tree':purity},
      'validation_contract':'After complete unchanged output clear, only localinvocation0 stores retainedStatus=0, then existing Workgroup-exec/Device-memory/2376 initial AllMemory join. Native loop initializes i=localinvocation, tests unsigned i<46, calls unchanged ReadOriginal(i,row), tests returned member4==0, calls unchanged atomicOR1 Reject helper only on invalid. No loop return; increment64 and full original loop/switch structure. All64 lanes reach new Workgroup/Workgroup/AcquireRelease|WorkgroupMemory264 join; only then loadstatus!=0 gives uniformreturn. Valid indices0..45 each execute once; idle46..63 skip call but join. maxincrement109 no uintoverflow. Allshape prechecks beforeloop and chart/spin/lambda/h/latermath/program/root/status/publication suffix retained under complete explicit typed ID mapping.',
      'normalization':'Function identities=(actual OpName, same-name occurrence), recursive exact type/constant/import definitions, named globals with separately exact typed declarations, local result IDs paired by definition order in all69helpers and1141keptmainrows. Complete ordered instruction tokens/operands compared; moved Reject body not matched by global firstoccurrence. Semantic decorations are full multiplicity-preserving multiset; only entry interface declaration-list ordering normalized as an unordered set of the same8IDs. OpName debug order and assembler Bound change are not runtime ABI. Mapping stores every local pair and static-signature correspondence.',
      'comparison_sizes':{'baseline_SPV_bytes':identity(W/'originals/canonical-retained_transport.spv')['bytes'],'candidate_SPV_bytes':identity(paths['.spv'])['bytes'],'baseline_barriers':dict(Counter(re.findall(r'OpControlBarrier (.+)',a.text))),'candidate_barriers':dict(Counter(re.findall(r'OpControlBarrier (.+)',b.text)))},
      'decision':'EMITTED_VALIDATION_CONTRACT_PASS_RUNTIME_COST_UNPROVEN','runtime_admitted':False,'production_adopted':False,
      'limits':['Proof composes with reviewed source uniform-shape/single-row/immutable-input assumptions; not arbitrary malformed-program or whole-scientific formal equivalence.','No new storage/cache/32-lane/packing/numeric-policy change; other34SPVs/39mainarrays/7tables/twoowners not rebuilt here and require exact latercomposition.','ReadOriginal reachability check excludes global writes/atomics/barriers and allows only Functionstorage writes; new RejectAtomicOr is outside that purecalltree. No claim about driver instruction work/scalarization/broadcast/bank behavior.','No provider/test/compiler/disassembler/validator/projecthelper/Git/network execution by reviewer; this Python reader only reads saved files and exact recorded Linux births. Original8controls+sole22gateABBA stillrequired.']}
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'receipt':pin(target),'mapping':pin(idmap),'readable_delta':pin(diff)},indent=2))


if __name__=='__main__':main()
