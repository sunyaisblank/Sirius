"""Remove compound-opcode packing, retaining every tuple and arithmetic body."""
from pathlib import Path
import collections,hashlib,json,re,sys
R=Path.cwd();W=Path(__file__).resolve().parent;D=R/'bin/linux-gcc/native-cpu-direct-compound/ispc-source'
assert sys.argv[1:] in ([],['--write'])
def pin(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def function(source,name):
 m=re.search(r'^static (?:noinline )?[^\n]+\b'+name+r'\([^{}]*\)\n\{',source,re.M);assert m,name
 depth=1;end=m.end()
 while depth:depth+=(source[end]=='{')-(source[end]=='}');end+=1
 return source[m.start():end]
basis_path=R/'out/native-cpu-arithmetic-closure/ispc-derivation.json';basis=json.loads(basis_path.read_text())
inputs={str(p.relative_to(R)):pin(p) for p in (basis_path,Path(__file__).resolve())}
stages={};total=0;extended=0
for slug,old in basis['stages'].items():
 p=R/old['derived_path'];assert pin(p)==old['derived'];source=p.read_text();candidate=source
 decl='''
// Direct compound opcode stores use the same private bank and original masks.
extern "C" void sirius_native_dense_subtract8(uniform uint32_t * uniform scratch,
    uniform uint32_t destination, uniform uint32_t a, uniform uint32_t b, uniform uint32_t active);
'''
 marker='typedef int32 int32_t;\n';assert candidate.count(marker)==1;candidate=candidate.replace(marker,marker+decl,1)
 original=function(source,'EvaluateOneStatic_0');header=original[:original.index('{')]
 names=[v.strip().split()[-1] for v in header[header.index('(')+1:header.rindex(')')].split(',')]
 dst,a,b,c,row,base,ctx=names;regs=re.search(r'uniform const uint32_t (\w+) = StaticRegisters;',original)[1]
 guard=re.search(r'    if \(StaticOp == 2U \|\| StaticOp == 4U\) \{.*?        return true;\n    \}',original,re.S)[0]
 new=f'''    if (StaticOp == 2U || StaticOp == 3U || StaticOp == 4U || StaticOp == 5U) {{
        if ({a} >= {regs} || {b} >= {regs})
            return false;
        if (StaticOp == 2U)
            RTAddStore_0({dst}, {a}, {b}, {base}, {ctx});
        else if (StaticOp == 3U)
            RTSubtractStore_0({dst}, {a}, {b}, {base}, {ctx});
        else if (StaticOp == 4U)
            RTMultiplyStore_0({dst}, {a}, {b}, {base}, {ctx});
        else
            RTDivideStore_0({dst}, {a}, {b}, {base}, {ctx});
        return true;
    }}
    if (StaticOp == 6U) {{
        if ({a} >= {regs})
            return false;
        RTSqrtStore_0({dst}, {a}, {a}, {base}, {ctx});
        return true;
    }}'''
 assert candidate.count(guard)==1;candidate=candidate.replace(guard,new,1)
 anchor=function(candidate,'RTMultiplyStore_0');helpers=''
 for name,target in (('RTSubtractStore_0','subtract'),('RTDivideStore_0','divide'),('RTSqrtStore_0','sqrt')):
  helpers+=f'''

// Same bank/source snapshots; no temporary bank, returned tuple or repacking.
static noinline void {name}(uniform uint32_t SiriusDestination,
    uniform uint32_t SiriusA, uniform uint32_t SiriusB, uniform uint32_t SiriusRowBase,
    varying KernelContext_0 * uniform SiriusContext)
{{
    sirius_native_dense_{target}8((uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],
        SiriusDestination, SiriusA, SiriusB, (uniform uint32_t)lanemask());
}}'''
 assert candidate.count(anchor)==1;candidate=candidate.replace(anchor,anchor+helpers,1)
 inverse=candidate.replace(helpers,'',1).replace(new,guard,1).replace(decl,'',1);assert inverse==source
 calls=lambda text:[line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
 assert calls(candidate)==calls(source)
 counts=collections.Counter(map(int,re.findall(r'if\(!EvaluateOneStatic_0<(\d+)U,',source)))
 total+=len(calls(source));extended+=sum(counts[k] for k in (3,5,6))
 target=D/(slug+'.ispc');data=candidate.encode()
 if target.exists():assert target.read_bytes()==data
 else:assert sys.argv[1:]==['--write'];target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
 inputs[str(p.relative_to(R))]=pin(p)
 stages[slug]={'basis_path':str(p.relative_to(R)),'basis':pin(p),'derived_path':str(target.relative_to(R)),'derived':pin(target),'opcode_sites':{str(k):v for k,v in sorted(counts.items())},'literal_instruction_calls_preserved':len(calls(source)),'inverse_restores_exact_AG_source':True,'changed_function':'EvaluateOneStatic_0','new_direct_stores':['RTSubtractStore_0','RTDivideStore_0','RTSqrtStore_0']}
assert total==9380 and extended==936
result={'stages':stages,'instructions':total,'additional_direct_opcode_sites':extended,'source_inputs':inputs,
 'scope':'AI extends the existing direct-bank opcode2/4 path to opcode3/5/6. Only three new store trampolines, the subtract ABI declaration and the static opcode dispatch branch change. All9380ordered literal tuples, original generic evaluator, original arithmetic/value helpers and nonstatic users remain byte-identical. Binary destination/a/b and unary destination/a guards have original outcomes; numerical invalid values are still stored and return true, as before, not turned into structural rejection. Complete operand snapshots precede writes, including reused destination aliases; current lanemask and private register layout remain. No rational precomputation, host-rounded constants, FP/ISA flags, precision, eligibility, requests, fields, trials, quartics, scheduler, queue or governor change. Fresh focused/full38/C8/affected original control5 and original frame remain required.'}
p=W/'ispc-derivation.json'
if p.exists():assert json.loads(p.read_text())==result
else:assert sys.argv[1:]==['--write'];p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'literal_instruction_calls_preserved':total,'new_direct_opcode_sites':extended,'compiled':False}))
