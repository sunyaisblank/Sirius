import hashlib, importlib.util, json, pathlib, struct, subprocess
root=pathlib.Path.cwd()
spec=importlib.util.spec_from_file_location('retained_program',root/'src/sirius/kernels/retained_program.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
original=m.compile_program

def decode(p):
 regs={}; roots=p['outputs']; operations=p['operations']; offsets=p['layer_offsets']; prefix=p['prefix_instructions']
 assert offsets[0]==0 and offsets[-1]==p['instructions'] and prefix in offsets
 assert p['registers']==612 and len(roots)==100
 projected=[];mixed=0
 for start,end in zip(offsets,offsets[1:]):
  assert 1<=end-start<=64
  layer=[operations[5*i:5*i+5] for i in range(start,end)]
  inputs=set();destinations=set();writes={}
  mixed+=len({op for op,*_ in layer})>1
  for i,(op,dst,a,b,c) in zip(range(start,end),layer):
   assert dst not in destinations and dst<612;destinations.add(dst)
   if op in [0,1]:
    args=[a,b,c]
    if op==1 and 45<=a<=52:
     assert i>=prefix;projected.append(i)
   else:
    dependencies=[a]+([] if op in [6,7,8] else [b])+([c] if op in [10,11] else [])
    assert all(d in regs for d in dependencies);inputs.update(dependencies)
    args=[regs[d] for d in dependencies]
   writes[dst]=hashlib.sha256(json.dumps([op,args],separators=(',',':')).encode()).hexdigest()
  assert inputs.isdisjoint(destinations)
  regs.update(writes)
 assert projected
 return {'root_hashes':[regs[r] for r in roots], 'instructions':p['instructions'],'registers':p['registers'],'outputs':len(roots),'prefix':prefix,'layers':len(offsets)-1,'mixed_layers':mixed,'projected_input_loads':projected,'words':3+len(roots)+len(operations)+len(offsets)}

reports=[]
for builder,prefix in [('build_endpoint_program',250),('build_schwarzschild_endpoint_program',202)]:
 m.compile_program=original; baseline=getattr(m,builder)(True)
 def homogeneous(outputs,parallel=False,prefix_outputs=None,homogeneous=False):
  return original(outputs,parallel,prefix_outputs,homogeneous=parallel)
 m.compile_program=homogeneous; trial=getattr(m,builder)(True)
 a,b=decode(baseline),decode(trial)
 assert a['root_hashes']==b['root_hashes'] and a['prefix']==b['prefix']==prefix
 assert a['instructions']==b['instructions'] and b['mixed_layers']==0
 reports.append({'builder':builder,'all_100_output_dag_hashes_identical':True,'baseline':a,'trial':b})
receipt={'source_baseline':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'scope':'Independent decoding of emitted complete layers, ordered operand DAGs and projected-input joins. Diagnostic in-memory one-argument trial only; no numerical/device/qualification claim.','reports':reports}
path=root/'out/endpoint-schedule-review/dag-identity.json';path.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps([{'builder':r['builder'],'same_100_roots':True,'baseline_layers':r['baseline']['layers'],'trial_layers':r['trial']['layers'],'added_words':r['trial']['words']-r['baseline']['words']} for r in reports]))
