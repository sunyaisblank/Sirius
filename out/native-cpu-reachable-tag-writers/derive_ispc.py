"""Move all twelve private tag updates before their original terminal returns."""
from pathlib import Path
import hashlib,json,re,sys
R=Path.cwd();W=Path(__file__).resolve().parent;D=R/'bin/linux-gcc/native-cpu-reachable-tag-writers/ispc-source'
assert sys.argv[1:] in ([],['--write'])
def pin(p):return dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def function(s,name):
 m=re.search(r'^static [^\n]+\b'+name+r'\(',s,re.M);assert m,name
 start=s.index('{',m.end());end=start+1;depth=1
 while depth:depth+=(s[end]=='{')-(s[end]=='}');end+=1
 return s[m.start():end]
def publish(p,data):
 if p.exists():assert p.read_bytes()==data,p
 else:assert sys.argv[1:]==['--write'];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
base=R/'out/native-cpu-normalization-tags/ispc-derivation.json';old=json.loads(base.read_text());inputs={str(base.relative_to(R)):pin(base),str(Path(__file__).resolve().relative_to(R)):pin(Path(__file__).resolve())}
stages={};total=0
for slug,entry in old['stages'].items():
 p=R/entry['derived_path'];assert pin(p)==entry['derived'];original=p.read_text();candidate=original;changes={}
 for name in ['ScratchWrite_0','ScratchWriteUniform_0','WriteValue_0','WriteValueUniform_0']:
  before=function(original,name);assert before.count('    return;\n')==1
  head,update=before.split('    return;\n');assert update.endswith('}') and '->SiriusNormalized_0[' in update
  assert 'return' not in head and update.count('->SiriusNormalized_0[')==1
  after=head+update[:-1]+'    return;\n}'
  assert after.index('->SiriusNormalized_0[')<after.index('    return;\n')
  assert candidate.count(before)==1;candidate=candidate.replace(before,after,1)
  changes[name]=dict(old_body_sha256=hashlib.sha256(before.encode()).hexdigest(),new_body_sha256=hashlib.sha256(after.encode()).hexdigest(),old_update_unreachable=True,new_update_before_only_return=True)
 inverse=candidate
 for name in changes:inverse=inverse.replace(function(candidate,name),function(original,name),1)
 assert inverse==original
 calls=lambda s:[line for line in s.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
 assert calls(candidate)==calls(original);total+=len(calls(original))
 dest=D/(slug+'.ispc');publish(dest,candidate.encode());inputs[str(p.relative_to(R))]=pin(p)
 stages[slug]=dict(basis_path=str(p.relative_to(R)),basis=pin(p),derived_path=str(dest.relative_to(R)),derived=pin(dest),writers=changes,literal_instruction_calls_preserved=len(calls(original)),whole_source_inverse_restores_AJ=True)
assert total==9380
result=dict(instructions=total,stages=stages,source_inputs=inputs,
 defect='AJ appended each of four tag writes after the generated unconditional return, in all three stages. Partial invalidation and complete refresh were unreachable. Earlier standalone predicate and arithmetic fixtures did not call these writers and did not establish coherence. AK reused the same raw3 and inherits this defect.',
 scope='Only move twelve private metadata writes before their original terminal returns. All five payload stores, masks/indices, canonical9380tuples, generic/runtime guards, original arithmetic, tags/predicate, complete initialization, public layouts and nonempty refinement behavior remain. Reuse exact AK C++ helper/adapter and old scalar3/pure admission. Focused actual-writer old/new regression and current whole-program controls required; no full scientific/performance/release claim.')
publish(W/'ispc-derivation.json',(json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(dict(reachable_writer_updates=12,unchanged_literal_calls=total,compiled=False)))
