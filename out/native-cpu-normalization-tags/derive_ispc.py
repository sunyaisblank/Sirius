"""Maintain one exact predicate per private register, with complete invalidation."""
from pathlib import Path
import hashlib,json,re,sys
R=Path.cwd();W=Path(__file__).resolve().parent;D=R/'bin/linux-gcc/native-cpu-normalization-tags/ispc-source'
assert sys.argv[1:] in ([],['--write'])
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def function(s,name):
 m=re.search(r'^static [^\n]+\b'+name+r'\(',s,re.M);assert m,name
 start=s.index('{',m.end());end=start+1;depth=1
 while depth:depth+=(s[end]=='{')-(s[end]=='}');end+=1
 return s[m.start():end]
def publish(p,data):
 if p.exists():assert p.read_bytes()==data,p
 else:assert sys.argv[1:]==['--write'];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
basis=R/'out/native-cpu-direct-compound/ispc-derivation.json';old=json.loads(basis.read_text())
inputs={str(p.relative_to(R)):pin(p) for p in (basis,Path(__file__).resolve(),W/'normalization_tag.isph')};stages={};count=0
for slug,item in old['stages'].items():
 p=R/item['derived_path'];assert pin(p)==item['derived'];original=p.read_text();candidate=original;changes=[]
 def replace(a,b):
  global candidate
  assert candidate.count(a)==1,a[:120];candidate=candidate.replace(a,b,1);changes.append((a,b))
 marker='typedef int32 int32_t;\n'
 decl='\n#include "'+str(W/'normalization_tag.isph')+'"\n'
 for op in ('add','multiply','subtract','divide','sqrt'):
  decl+='extern "C" void sirius_native_tagged_'+op+'8(uniform uint32_t * uniform scratch,\n    uniform uint32_t * uniform tags, uniform uint32_t destination, uniform uint32_t a,\n    uniform uint32_t b, uniform uint32_t active);\n'
 replace(marker,marker+decl)
 field,=re.findall(r'    uint32_t retainedScratch_0\[(\d+)\];',candidate);capacity=int(field)//5;assert int(field)==5*capacity
 before='    uint32_t retainedScratch_0['+field+'];';replace(before,before+'\n    uint32_t SiriusNormalized_0['+str(capacity)+'];')
 writers={}
 for name in ('ScratchWrite_0','ScratchWriteUniform_0'):
  before=function(candidate,name);header=before[:before.index('{')];args=[x.strip().split()[-1] for x in header[header.index('(')+1:header.rindex(')')].split(',')];index=args[0];ctx=args[-1]
  assert before.count('->retainedScratch_0[')==1
  after=before[:-1]+'    '+ctx+'->SiriusNormalized_0['+index+'] = 0U;\n}'
  replace(before,after);writers[name]='Every partial payload write invalidates this same register tag under the original active mask.'
 for name in ('WriteValue_0','WriteValueUniform_0'):
  before=function(candidate,name);header=before[:before.index('{')];args=[x.strip().split()[-1] for x in header[header.index('(')+1:header.rindex(')')].split(',')];index,value,_,ctx=args
  assert len(re.findall(r'ScratchWrite(?:Uniform)?_0\(',before))==5
  update=f'''    {ctx}->SiriusNormalized_0[{index}] = SiriusNormalizationTag_0(
        RPBits_0({value}->hi_0), RPBits_0({value}->lo_0), RPBits_0({value}->tail_0),
        RPBits_0({value}->error_0), {value}->valid_0);
'''
  replace(before,before[:-1]+update+'}');writers[name]='All five original payload words are written first, then the exact shared integer predicate is computed from that same complete tuple.'
 for op,stem in (('add','RTAddStore_0'),('multiply','RTMultiplyStore_0'),('subtract','RTSubtractStore_0'),('divide','RTDivideStore_0'),('sqrt','RTSqrtStore_0')):
  before=function(candidate,stem)
  a='sirius_native_dense_'+op+'8((uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],'
  b='sirius_native_tagged_'+op+'8((uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],\n        (uniform uint32_t * uniform)&SiriusContext->SiriusNormalized_0[0],'
  assert before.count(a)==1;replace(before,before.replace(a,b,1));writers[stem]='The new CPP entry snapshots both operands/tags, uses the original fast domain only with exact current tags, and refreshes the active destination tag from its actual result.'
 # No payload writes bypass the two audited low-level accessors, apart from
 # these five C++ direct stores. Context initialization includes the new field.
 assert len(re.findall(r'->retainedScratch_0\[[^\n]*\]\s*=',original))==2
 assert 'memset((void * uniform)&context, 0, sizeof(context));' in original
 assert len(re.findall(r'->SiriusNormalized_0\[[^\n]*\] = 0U;',candidate))==2
 literal=lambda s:[line for line in s.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
 assert literal(candidate)==literal(original);count+=len(literal(original))
 inverse=candidate
 for a,b in reversed(changes):assert inverse.count(b)==1;inverse=inverse.replace(b,a,1)
 assert inverse==original
 destination=D/(slug+'.ispc');publish(destination,candidate.encode());inputs[str(p.relative_to(R))]=pin(p)
 stages[slug]={'basis_path':str(p.relative_to(R)),'basis':pin(p),'derived_path':str(destination.relative_to(R)),'derived':pin(destination),'literal_instruction_calls_preserved':len(literal(original)),'register_tag_capacity':capacity,'added_private_metadata_bytes_per_packet':capacity*8*4,'writer_protocol':writers,'whole_source_inverse_restores_AI':True}
assert count==9380
result={'instructions':count,'stages':stages,'source_inputs':inputs,
 'scope':'AJ caches only the unchanged normalized-input predicate alongside the same complete private five-word register. The metadata is initialized with the whole context. Both generic/uniform low-level writes invalidate it, including validity-only clears; both complete tuple writes recompute it after all payload words. Five direct CPP stores also refresh it after writes. All literal tuples, generic evaluation/guards, math/value helpers and public input/output layouts remain. No assumed normalization theorem, literal/physical arithmetic folding, extra limb, persistent/input cache, shader-cache policy change, stronger eligibility, altered precision/FP/ISA/governors/requests/fields/RHS/trials/quartics/order/refusal/rollback. Source derivation only; shared ISPC predicate fixture and focused arithmetic/tag/alias checks, actual code, full38/C8/affected control5/frame required.'}
publish(W/'ispc-derivation.json',(json.dumps(result,indent=2)+'\n').encode())
print(json.dumps({'original_literal_calls':count,'stages':3,'all_writes_invalidate_or_refresh':True,'compiled':False}))
