"""Move only original AM identity stores ahead of its C ABI boundary."""
from pathlib import Path
import hashlib,json,re,sys
R=Path.cwd();W=Path(__file__).resolve().parent;D=R/'bin/linux-gcc/native-cpu-local-identities/ispc-source'
assert sys.argv[1:]==['--write']
def pin(p):return dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def publish(p,data):
 if p.exists():assert p.read_bytes()==data
 else:p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
oldpath=R/'out/native-cpu-value-fact-stores/ispc-derivation.json';old=json.loads(oldpath.read_text())
helper=W/'identity_store.isph';inputs={str(p.relative_to(R)):pin(p) for p in (oldpath,helper,Path(__file__).resolve(),R/'out/native-cpu-value-fact-stores/retained8.cpp')};stages={};total=0
for slug,row in old['stages'].items():
 p=R/row['derived_path'];assert pin(p)==row['derived'];original=p.read_text();candidate=original
 anchor='#include "'+str(R/'out/native-cpu-value-fact-stores/value_facts.isph')+'"'
 include='#include "'+str(helper)+'"';assert candidate.count(anchor)==1
 candidate=candidate.replace(anchor,anchor+'\n'+include,1);inserts=[]
 for name,op in [('Add',2),('Subtract',3),('Multiply',4)]:
  match=re.search(r'^static noinline void RT'+name+r'Store_0\([^}]+\n\}',candidate,re.M);assert match
  body=match[0];a=body.index('{')+1
  inserted='\n    if(SiriusTryIdentity_0((uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],\n        (uniform uint32_t * uniform)&SiriusContext->SiriusNormalized_0[0],\n        SiriusDestination,SiriusA,SiriusB,'+str(op)+'U)) return;'
  changed=body[:a]+inserted+body[a:];candidate=candidate.replace(body,changed,1);inserts.append(inserted)
 inverse=candidate.replace('\n'+include,'',1)
 for insert in inserts:assert inverse.count(insert)==1;inverse=inverse.replace(insert,'',1)
 assert inverse==original
 calls=lambda s:[x for x in s.splitlines() if x.strip().startswith('if(!EvaluateOneStatic_0<')]
 assert calls(candidate)==calls(original);total+=len(calls(original));dest=D/(slug+'.ispc');publish(dest,candidate.encode())
 stages[slug]=dict(basis_path=str(p.relative_to(R)),basis=pin(p),derived_path=str(dest.relative_to(R)),derived=pin(dest),literal_instruction_calls_preserved=len(calls(original)),whole_source_inverse_restores_AM=True,unchanged_fallback_wrapper_bodies=True)
 inputs[str(p.relative_to(R))]=pin(p)
assert total==9380
proof=dict(instructions=total,stages=stages,source_inputs=inputs,scope='One shared integer-only original AM identity/refusal store before the same three Add/Subtract/Multiply C ABI calls in all three stages. All active lanes must have complete coherent metadata and select an original identity/refusal before any store. All payload snapshots precede writes, zero/refusal-only packets load no payload, and unknown/legacy/general packets enter exact unchanged AM helper. No new array/context allocation, arithmetic expression, table/tuple/layer/mask/bounds/input/output/FP/ISA/count/cache/governor/scheduler change; all writers and their tested reachable fact protocol are byte identical. Full inverse restores AM. Actual fixture/code/full38/C8 and original adoption gates remain required.')
publish(W/'ispc-derivation.json',(json.dumps(proof,indent=2)+'\n').encode())
fixture='typedef uint32 uint32_t;\n#include "'+str(helper)+'"\n'
for name,op in [('add',2),('subtract',3),('multiply',4)]:
 fixture+='''extern "C" void sirius_native_tagged_NAME8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask);
export void local_tagged_NAME8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask) {
    if((mask&(1U<<programIndex))!=0U) {
        if(!SiriusTryIdentity_0(bank,tags,d,a,b,OPU))
            sirius_native_tagged_NAME8(bank,tags,d,a,b,(uniform uint32_t)lanemask());
    }
}
'''.replace('NAME',name).replace('OP',str(op))
publish(W/'fixture.ispc',fixture.encode());print(json.dumps(dict(derived_stages=3,unchanged_literal_calls=total,compiled=False)))
