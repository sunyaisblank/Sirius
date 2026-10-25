"""Extract actual old/new writer bodies and exact stage context declarations."""
from pathlib import Path
import hashlib,json,re,runpy,sys
R=Path.cwd();W=Path(__file__).resolve().parent
assert sys.argv[1:]==['--write']
namespace=runpy.run_path(str(W/'derive_ispc.py'),run_name='writer_derivation')
function=namespace['function'];pin=namespace['pin'];J=lambda p:json.loads(p.read_text())
derived=J(W/'ispc-derivation.json');pred=W/'value_facts.isph'
source='typedef uint32 uint32_t;\n#define F32_asuint(x) intbits(x)\n'+pred.read_text()+'\nstruct GlobalParams_0 { uniform uint32_t unused; };\n'
triple=None;proofs=[]
for slug,stage in derived['stages'].items():
 for label,key in [('old','counterexample_basis_path'),('fixed','derived_path')]:
  text=(R/stage[key]).read_text();prefix=slug+'_'+label+'_'
  declaration=re.search(r'struct RetainedTriple_0\n\{.*?\n\};',text,re.S)[0]
  if triple is None:triple=declaration;source+=triple+'\n'
  else:assert declaration==triple
  context=re.search(r'struct KernelContext_0\n\{.*?\n\};',text,re.S)[0]
  names=['RPBits_0','ScratchAddress_0','ScratchAddressUniform_0','ScratchWrite_0','ScratchWriteUniform_0','WriteValue_0','WriteValueUniform_0']
  bodies={name:function(text,name) for name in names}
  renames=['KernelContext_0',*names]
  def rename(s):
   for name in renames:s=re.sub(r'\b'+name+r'\b',prefix+name,s)
   return s
  source+=rename(context)+'\n'
  for name,body in bodies.items():
   changed=rename(body)
   inverse=changed
   for old in renames:inverse=re.sub(r'\b'+prefix+old+r'\b',old,inverse)
   assert inverse==body;source+=changed+'\n'
  source+=f'''export void writer_{slug}_{label}(uniform uint32_t* uniform bank,
    uniform uint32_t* uniform tags, uniform uint32_t* uniform value,
    uniform uint32_t operation, uniform uint32_t index, uniform uint32_t component,
    uniform uint32_t mask, uniform uint32_t vary_index) {{
    {prefix}KernelContext_0 context;
    for(uniform uint32_t r=0;r<3U;++r) {{
        for(uniform uint32_t c=0;c<5U;++c)
            context.retainedScratch_0[r*5U+c]=bank[(r*5U+c)*programCount+programIndex];
        context.SiriusNormalized_0[r]=tags[r*programCount+programIndex];
    }}
    if((mask&(1U<<programIndex))!=0U) {{
        uint32_t varying_index=vary_index!=0U ? (index+programIndex)%3U : index;
        RetainedTriple_0 v;
        v.hi_0=floatbits(value[programIndex]);v.lo_0=floatbits(value[programCount+programIndex]);
        v.tail_0=floatbits(value[2U*programCount+programIndex]);v.error_0=floatbits(value[3U*programCount+programIndex]);
        v.valid_0=value[4U*programCount+programIndex];
        if(operation==0U) {prefix}ScratchWrite_0(varying_index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==1U) {prefix}ScratchWriteUniform_0(index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==2U) {prefix}WriteValue_0(varying_index,&v,0U,&context);
        else {prefix}WriteValueUniform_0(index,&v,0U,&context);
    }}
    for(uniform uint32_t r=0;r<3U;++r) {{
        for(uniform uint32_t c=0;c<5U;++c)
            bank[(r*5U+c)*programCount+programIndex]=context.retainedScratch_0[r*5U+c];
        tags[r*programCount+programIndex]=context.SiriusNormalized_0[r];
    }}
}}
'''
  proofs.append(dict(stage=slug,version=label,stage_source=stage[key],stage_identity=pin(R/stage[key]),exact_context_declaration=hashlib.sha256(context.encode()).hexdigest(),writer_body_SHA256={n:hashlib.sha256(bodies[n].encode()).hexdigest() for n in names},source_inverse_alpha_renames_exact=True))
path=W/'writers.ispc';assert not path.exists();path.write_text(source)
proof=dict(fixture=pin(path),extracts=proofs,predicate=dict(path=str(pred.relative_to(R)),**pin(pred)),scope='Six exports call exact alpha-renamed old/fixed writer/address/bitcast bodies in exact stage KernelContext declarations. GlobalParams is an unused opaque pointer target; only first three initialized register/tag slots are accessed. Every active/inactive word in those slots is returned. No production arithmetic copied or modeled. Original scalar and current complete-fact arithmetic supply the independent stale-tag integration witness; actual raw provider/code/source joins remain separate.')
(W/'writer-fixture-source.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(dict(extracted_versions=6,writer_bodies=24,compiled=False)))
