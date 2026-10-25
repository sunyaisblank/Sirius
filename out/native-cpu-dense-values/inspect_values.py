"""Require unchanged prior math/adapter code and both value ABI references."""
from pathlib import Path
import json
import subprocess

R=Path.cwd()
W=Path(__file__).resolve().parent
doc=lambda p:json.loads(p.read_text())
compiled=doc(W/'ispc-values-compile-result.json')
owner=doc(W/'compile/ispc-values/owner.json')
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert compiled['all_objects_compiled'] and compiled['new_compilations']==3
assert compiled['new_CPP_compilations']==0 and compiled['new_relocatable_links']==1
proof=compiled['adapter_input_code_association']
assert proof['passed'] and all(v['all_retained_function_opcodes_match'] for v in proof['objects'])
assert doc(W/'object-code-inspection-values.json')['objects'] and doc(W/'bridge-code-inspection.json')['passed']
for row in compiled['objects']:
    relocations=subprocess.check_output(['/usr/bin/objdump','-r',str(R/row['object']['path'])],text=True)
    for name in ('sirius_native_dense_values_add8','sirius_native_dense_values_multiply8'):assert name in relocations
result={'passed':True,'prior_math_and_new_packing_code_association':proof,
        'both_value_ABI_references_in_all3':True,
        'scope':'Actual object/ISA and unchanged prior math/adapter opcode associations. New packing-only functions add no arithmetic source or opcode changes to that prior object. Actual relocated targets/full38/C8/frame still required; no renderer/scientific qualification.'}
with (W/'values-code-inspection.json').open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps({'passed':True,'prior_dense_math_code_unchanged':True,'stages':3}))
