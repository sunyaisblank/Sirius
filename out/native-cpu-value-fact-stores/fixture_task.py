"""Compile/link or execute the sealed focused writer regression exactly once."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'fixture';J=lambda p:json.loads(p.read_text())
def pin(p):
 with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())
mode,=sys.argv[1:];p=J(F/'prepared.json');exe=R/p['executable']
def verify():
 for name,v in p['inputs'].items():assert pin(R/name)==v,name
verify()
if mode=='build':
 assert not exe.parent.exists() and not (F/'build-result.json').exists();exe.parent.mkdir(parents=True)
 r=dict(passed=False,commands=[],new_ISPC_compilations=0,new_CPP_compilations=0,new_links=0)
 def save():(F/'build-result.json').write_text(json.dumps(r,indent=2)+'\n')
 save()
 for i,args in enumerate(p['commands']):
  row=dict(argv=args,started=True);r['commands'].append(row);save();code=subprocess.run(args,cwd=R).returncode;row['exit']=code
  r[['new_ISPC_compilations','new_CPP_compilations','new_links'][i]]+=1;save()
  if code:raise SystemExit(code)
  output=Path(args[-1]);row['product']=dict(path=str(output.relative_to(R)),**pin(output));save()
 sys.path.insert(0,str(W));from elf_writer_proof import associate_fixture
 objects=[r['commands'][0]['product'],r['commands'][1]['product'],p['helper']]
 proof=associate_fixture(exe,objects);assert proof['passed']
 (F/'linked-code.json').write_text(json.dumps(proof,indent=2)+'\n')
 text=subprocess.check_output(['ldd','-r',str(exe)],text=True,stderr=subprocess.STDOUT);assert 'undefined symbol' not in text and 'not found' not in text
 r.update(executable=dict(path=p['executable'],**pin(exe)),objects=objects,linked_code=pin(F/'linked-code.json'),
  runtime_libraries={str(Path(n).resolve()):pin(Path(n).resolve()) for n in re.findall(r'(/[^\s()]+)',text) if Path(n).is_file()})
 verify();r['passed']=True;save();print(json.dumps(dict(build_passed=True,compiled_fixture_only=True)))
else:
 assert mode=='run' and not (F/'result.json').exists();built=J(F/'build-result.json');assert built['passed'] and pin(exe)=={k:built['executable'][k] for k in ('bytes','sha256')}
 output=subprocess.check_output([str(exe)],text=True);r=json.loads(output)
 assert r['passed'] and r['old_stale_tag_arithmetic_counterexamples']==24 and r['case_tuples']==49 and r['stages']==3 and r['writer_forms_per_stage']==4
 verify();r.update(executable=built['executable'],linked_code=pin(F/'linked-code.json'),scope='Actual copied old/fixed writer bodies only; 24legacy/rich stale-tag arithmetic counterexamples and complete corrected word/tag checks. Source/code/provider and full renderer qualification remain separately required.')
 (F/'result.json').write_text(json.dumps(r,indent=2)+'\n');print(output,end='')
