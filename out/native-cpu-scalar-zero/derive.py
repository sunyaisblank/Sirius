"""Put the reviewed exact-zero arithmetic in the actually used D scalar path.

No field packing or row scheduling changes. New full scalar source/header
copies use the original exports and exactly the original compiler recipe.
"""
from pathlib import Path
import difflib
import hashlib
import json
import re
import sys

R=Path.cwd()
W=Path(__file__).resolve().parent
B=R/'bin/linux-gcc/native-cpu-scalar-zero/typed-batch'
assert sys.argv[1:] in ([],['--write'])


def pin(data):
    return {'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}


def publish(path,data):
    if path.exists():
        assert path.read_bytes()==data,path
    else:
        assert sys.argv[1:]==['--write'],path
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)


base=json.loads((R/'out/native-cpu-static-program/header-preparation.json').read_text())
zero=json.loads((R/'out/native-cpu-zero-term/source-proof.json').read_text())
record={'instructions':0,'stages':{},'source_inputs':{},
    'scope':'Source derivation only. Apply the previously derived exact signed-zero arm and explicit Term inlining to D scalar arithmetic. Every other original header/wrapper byte, all9380instructions, scalar row control, coefficients, output publication and fallback ABI remain unchanged. No compile/runtime/benefit/scientific qualification.',
    'premise':'Q/R/S small Transport grids select D scalar rows, which never received the R/S ISPC-only Term changes. Field-SIMD T first original candidate benefit slot timed out at300s; no frame or unchanged retry is justified. This narrower scalar candidate needs its own actual build/code/complete fixture/C8/ABBA/all8/frame acceptance.',
    'semantic_argument':zero['semantic_argument'],'signed_zero_table':zero['signed_zero_table']}
for path in (Path(__file__).resolve(),R/'out/native-cpu-static-program/header-preparation.json',
             R/'out/native-cpu-zero-term/source-proof.json'):
    record['source_inputs'][str(path.relative_to(R))]=pin(path.read_bytes())
patch=''
for slug in ('transport','endpoint','dopri_phase'):
    original=R/base['headers'][slug]['derived_path']
    text=original.read_text()
    assert pin(text.encode())==base['headers'][slug]['derived']
    old='static RetainedTriple_0 RTAccumulateTerm_0('
    new='static inline __attribute__((always_inline)) RetainedTriple_0 RTAccumulateTerm_0('
    assert text.count(old)==1
    start=text.index(old)
    refusal='        return RTInvalid_0();\n    }'
    end=text.index(refusal,start)+len(refusal)
    arm=zero['stages'][slug]['zero_arm']
    pointer=re.search(r'RTAccumulateTerm_0\(RetainedTriple_0 \* (\w+), float (\w+)\)',text)
    assert pointer and pointer[1]+'->hi_0' in arm and 'RPBits_0('+pointer[2]+')' in arm
    candidate=(text[:end]+arm+text[end:]).replace(old,new,1)
    assert candidate.replace(new,old,1).replace(arm,'',1)==text
    literal=lambda value:[line for line in value.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert literal(candidate)==literal(text)
    instructions=len(literal(text)); record['instructions']+=instructions
    header=B/(slug+'.native.hpp')
    wrapper=R/('bin/linux-gcc/native-cpu-slang-bridge/typed-batch/'+slug+'.cpp')
    source=B/(slug+'.cpp')
    publish(header,candidate.encode());publish(source,wrapper.read_bytes())
    patch+=''.join(difflib.unified_diff(text.splitlines(True),candidate.splitlines(True),
        fromfile='a/'+str(original.relative_to(R)),tofile='b/'+str(header.relative_to(R))))
    record['stages'][slug]={'basis_header':{'path':str(original.relative_to(R)),**pin(text.encode())},
        'header':{'path':str(header.relative_to(R)),**pin(candidate.encode())},
        'source':{'path':str(source.relative_to(R)),**pin(wrapper.read_bytes())},
        'literal_instruction_calls_retained':instructions,
        'inverse_removes_only_zero_arm_and_inline_attribute':True,
        'nonzero_body_unchanged':True,'wrapper_bytes_identical':True}
    for path in (original,wrapper):
        record['source_inputs'][str(path.relative_to(R))]=pin(path.read_bytes())
assert record['instructions']==9380
publish(W/'scalar-zero.patch',patch.encode())
publish(W/'derivation.json',(json.dumps(record,indent=2)+'\n').encode())
print(json.dumps({'new_scalar_sources':3,'new_scalar_headers':3,'instructions':9380,'compiled':False}))
