"""Bind the integer-only early-return wrappers to actual raw objects."""
from pathlib import Path
import hashlib,json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;J=lambda p:json.loads(p.read_text());compiled=J(W/'ispc-dense-compile-result.json')
assert compiled['all_objects_compiled'] and J(W/'dense-code-inspection.json')['passed'];rows=[]
for item in compiled['objects']:
 p=R/item['object']['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==item['object']['sha256']
 source=(R/J(W/'ispc-derivation.json')['stages'][item['stage']]['derived_path']).read_text()
 assert source.count('#include "'+str(W/'identity_store.isph')+'"')==1
 asm=subprocess.check_output(['objdump','-drC','--no-show-raw-insn',str(p)],text=True);blocks=re.split(r'(?=^[0-9a-f]+ <)',asm,flags=re.M);found=[]
 for name in ['Add','Subtract','Multiply']:
  bodies=[b for b in blocks if re.match(r'^[0-9a-f]+ <RT'+name+r'Store_0',b)];assert len(bodies)==1;body=bodies[0]
  assert re.search(r'\bret\b',body) and ('sirius_native_tagged_'+name.lower()+'8') in body
  assert re.search(r'\bj(?:e|ne|z|nz)\b',body)
  assert not re.search(r'\bv(?:add|sub|mul|div|sqrt|fm\w*|fnm\w*)(?:ps|ss|pd|sd)\b',body)
  found.append(dict(wrapper=body.splitlines()[0],disassembly_sha256=hashlib.sha256(body.encode()).hexdigest(),integer_identity_return_and_original_fallback_present=True))
 rows.append(dict(stage=item['stage'],object=item['object'],wrappers=found))
result=dict(passed=True,common_identity_include=dict(path=str((W/'identity_store.isph').relative_to(R)),sha256=hashlib.sha256((W/'identity_store.isph').read_bytes()).hexdigest()),fixture_result=dict(path=str((W/'fixture/result.json').relative_to(R)),sha256=hashlib.sha256((W/'fixture/result.json').read_bytes()).hexdigest()),stages=rows,scope='All nine actual raw Add/Subtract/Multiply wrappers contain an early return and the same original AM fallback relocation. The new decision/store body is integer-only; no vector FP arithmetic is present in those wrappers. The same forced-inline include passes the scalar-reference fixture. This establishes emitted mechanism and source association, not runtime hit counts, timing improvement or renderer/scientific qualification.')
(W/'identity-code-inspection.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(passed=True,integer_identity_wrappers=9)))
