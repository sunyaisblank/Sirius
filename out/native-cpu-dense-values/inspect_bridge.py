"""Inspect the new packing-only object; arithmetic is in the unchanged AC object."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R=Path.cwd()
W=Path(__file__).resolve().parent
doc=lambda p:json.loads(p.read_text())
build=doc(W/'bridge-build-result.json')
owner=doc(W/'bridge/build/owner.json')
assert build['passed'] and owner['passed'] and owner['status']=='terminal'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
row=build['object'];path=R/row['path']
assert hashlib.file_digest(path.open('rb'),'sha256').hexdigest()==row['sha256']
undefined=subprocess.check_output(['/usr/bin/nm','-u',str(path)],text=True)
for name in ('sirius_native_dense_add8','sirius_native_dense_multiply8'):assert name in undefined
asm=subprocess.check_output(['/usr/bin/objdump','-d','--no-show-raw-insn',str(path)],text=True)
ops=re.findall(r'^\s*[0-9a-f]+:\s+(\w+)',asm,re.M)
forbidden={op for op in ops if op.startswith(('v','f','addp','subp','mulp','divp','sqrtp')) or
           op in ('ldmxcsr','stmxcsr','addss','subss','mulss','divss','sqrtss')}
assert not forbidden,forbidden
result={'passed':True,'object':row,'undefined_symbols':undefined.splitlines(),
        'mnemonics':sorted(set(ops)),'forbidden_vector_math_or_FP_control':sorted(forbidden),
        'scope':'New baseline packing object only: initialized buffer copies and calls to unchanged AC dense math, no arithmetic implementation or ISA expansion. Complete partial-link opcode joins and actual linked targets plus full38/C8/frame still required.'}
with (W/'bridge-code-inspection.json').open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps({'passed':True,'packing_only_CPP_object':True}))
