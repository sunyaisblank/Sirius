"""Restore original per-lane terminal returns when no refinement lane remains."""
from pathlib import Path
import hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent
base=R/'out/native-cpu-normalization-tags/retained8.cpp';original=base.read_text()
def pin(p):return dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
source=original
divide_return='    if(!any(work)) return choose(identity,a,invalid());\n'
sqrt_return='    if(!any(work)) return choose(input_invalid,a,choose(exact_zero,constant_zero(),invalid()));\n'
changes=[
 ('    U work=eligible&~identity;\n',divide_return),
 ('    Triple residual=subtract_triple(a,multiply_triple(b,q,work),work);\n',divide_return),
 ('    work=work&nz(residual.valid);\n    Triple result=pack(q.hi,q.lo,q.tail,ratio_upper(absolute_upper(residual),lower));\n',divide_return),
 ('           ~less_equal(absolute_lower(a),zero());\n',sqrt_return),
 ('    work=work&nz(q.valid)&~less_equal(q.hi,zero())&~less_equal(lower,zero());\n',sqrt_return),
 ('    work=work&nz(residual.valid);\n    Triple result=pack(q.hi,q.lo,q.tail,ratio_upper(absolute_upper(residual),lower));\n',sqrt_return),
]
for i,(anchor,addition) in enumerate(changes):
 if i==1:
  # Two identical residual expressions occur: loop body and final residual.
  # Select only the original four-space final expression (loop has eight).
  needle='\n'+anchor;assert source.count(needle)==1
  source=source.replace(needle,'\n'+addition+anchor,1)
 elif i in (2,5):
  # This anchor exists once per operation; restrict replacement to its body.
  start=source.index('inline Triple '+('divide_triple' if i==2 else 'sqrt_triple'))
  end=source.index('\n}\n',start)+3;body=source[start:end];assert body.count(anchor)==1
  body=body.replace(anchor,anchor.replace('    Triple result=',addition+'    Triple result='),1)
  source=source[:start]+body+source[end:]
 else:
  assert source.count(anchor)==1;source=source.replace(anchor,anchor+addition,1)
assert source.count(divide_return)==3 and source.count(sqrt_return)==3
restored=source.replace(divide_return,'').replace(sqrt_return,'');assert restored==original
target=W/'retained8.cpp';assert not target.exists();target.write_text(source)
proof=dict(basis=dict(path=str(base.relative_to(R)),**pin(base)),candidate=dict(path=str(target.relative_to(R)),**pin(target)),
 new_uniform_terminal_returns=6,exact_inverse_restores_AJ=True,unchanged_refinement_count=3,
 argument=[
 'Divide identity is the original eligible&exact-denominator-one mask. Original final selection returns a for those lanes and canonical invalid for every other lane when work is empty. Both masks only lose active lanes. The three inserted exits return precisely that same five-word selection before seed, final residual or final pack.',
 'Sqrt original final selection returns original a verbatim for invalid input, canonical valid zero for exact-zero input and canonical invalid for remaining lanes when work is empty. All three exits retain that original precedence, arbitrary invalid words and signed-zero behavior.',
 'No active refinement lane skips arithmetic: any nonempty work mask enters the identical original seed, three ordered correction steps, final residual, radius and normalization/refusal code. Mixed lanes retain the original vector computation and selection. Destination alias snapshots, masks, complete bank writes and normalization-tag refresh are unchanged.',
 'The provider already masks IEEE32 exceptions and restores the complete saved worker state. Fewer unused calculations may change intermediate exception status; no preservation of unobserved intermediate status is claimed. FP controls/ISA/rounding and the original scalar mathematical contract remain.',
 'This is source reasoning plus an exact byte inverse, not machine-checked equivalence, runtime benefit, full science or release acceptance. Existing 49-case complete-bank/tag/mask/alias fixture is reused, with broadcast identity/refusal/zero cases. Actual code joins and original whole-path checks remain required.'
 ],scope='Six empty-work terminal returns in existing C++ Divide/Sqrt only. Same five public words/private normalized tags, raw3 objects, canonical9380tuples/programs/SPVs, eligibility, refinement formulas/order, all fields/RHS/trials/quartics, scheduler, queue, counts, cache/governors and complete fallback.')
(W/'cpp-source-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(dict(source_derived=True,empty_work_returns=6,compiled=False)))
