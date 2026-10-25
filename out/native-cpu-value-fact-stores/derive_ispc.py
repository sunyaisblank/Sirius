"""Populate complete value facts in the corrected, reachable AL writers."""
from pathlib import Path
import hashlib,json,re,sys
R=Path.cwd();W=Path(__file__).resolve().parent;D=R/'bin/linux-gcc/native-cpu-value-fact-stores/ispc-source'
assert sys.argv[1:] in ([],['--write'])
def pin(p):return dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def function(s,name):
 m=re.search(r'^static [^\n]+\b'+name+r'\(',s,re.M);assert m,name
 a=s.index('{',m.end());b=a+1;depth=1
 while depth:depth+=(s[b]=='{')-(s[b]=='}');b+=1
 return s[m.start():b]
def publish(p,data):
 if p.exists():assert p.read_bytes()==data,p
 else:assert sys.argv[1:]==['--write'];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
base=R/'out/native-cpu-reachable-tag-writers/ispc-derivation.json';old=json.loads(base.read_text());inputs={str(p.relative_to(R)):pin(p) for p in [base,Path(__file__).resolve(),W/'value_facts.isph',R/'out/native-cpu-normalization-tags/normalization_tag.isph']};stages={};total=0
for slug,row in old['stages'].items():
 p=R/row['derived_path'];assert pin(p)==row['derived'];original=p.read_text()
 a='#include "'+str(R/'out/native-cpu-normalization-tags/normalization_tag.isph')+'"'
 b='#include "'+str(W/'value_facts.isph')+'"';assert original.count(a)==1
 candidate=original.replace(a,b,1);assert candidate.count('SiriusNormalizationTag_0(')==2
 candidate=candidate.replace('SiriusNormalizationTag_0(','SiriusValueFacts_0(')
 assert candidate.replace(b,a,1).replace('SiriusValueFacts_0(','SiriusNormalizationTag_0(')==original
 for name in ['ScratchWrite_0','ScratchWriteUniform_0','WriteValue_0','WriteValueUniform_0']:
  body=function(candidate,name);assert body.count('    return;\n')==1
  assert body.index('->SiriusNormalized_0[')<body.index('    return;\n')
  if name.startswith('Scratch'):assert re.search(r'->SiriusNormalized_0\[[^\n]+\] = 0U;',body)
  else:assert body.count('SiriusValueFacts_0(')==1
 calls=lambda s:[line for line in s.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
 assert calls(candidate)==calls(original);total+=len(calls(original))
 dest=D/(slug+'.ispc');publish(dest,candidate.encode());inputs[str(p.relative_to(R))]=pin(p)
 counterexample=R/row['basis_path'];assert pin(counterexample)==row['basis']
 stages[slug]=dict(basis_path=str(p.relative_to(R)),basis=pin(p),derived_path=str(dest.relative_to(R)),derived=pin(dest),counterexample_basis_path=str(counterexample.relative_to(R)),counterexample_basis=pin(counterexample),literal_instruction_calls_preserved=len(calls(original)),whole_source_inverse_restores_AL=True,all_four_writers_update_before_only_return=True)
assert total==9380
result=dict(instructions=total,stages=stages,source_inputs=inputs,
 scope='Replace only the shared metadata include and two complete-write predicate calls per stage. Existing one-word private metadata now carries exact normalization/identity/nonzero-valid facts and completeness; all partial writes still clear before their original returns. Every complete tuple refresh is reachable after five payload stores. Same arrays/offsets/published words, masks, canonical9380tuples, arithmetic, guards, precision/count/cache/governors. Exact inverse restores corrected AL. Original AJ bad writers retained only for counterexamples, never current runtime.')
publish(W/'ispc-derivation.json',(json.dumps(result,indent=2)+'\n').encode());print(json.dumps(dict(derived_stages=3,unchanged_literal_calls=total,compiled=False)))
