"""Inspect the retained Transport object; compile remaining Endpoint once.

The failed Transport receipt expected IR and native static call-site counts to
match. That assumption is stronger than the compiler contract: equivalent
control-flow tails may share one machine call without dropping logical work.
Keep that failed receipt. This continuation never recompiles Transport.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

R=Path.cwd();W=Path(__file__).resolve().parent;O=R/'bin/linux-gcc/native-cpu-register-value'
V=R/'out/native-cpu-strict-ir-linking';D=R/'bin/linux-gcc/native-cpu-strict-ir-linking/toolchain/bin'
PREFIX='sirius_register_value_'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())


def doc(path):return json.loads(path.read_text())


def write(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')


def native(stage,path):
    text=subprocess.check_output(['/usr/bin/objdump','-drw',str(path)],text=True)
    instructions=re.findall(r'^\s*[0-9a-f]+:\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9.]*)',text,re.M)
    assert instructions and not [x for x in instructions if x.startswith(('vfm','vfnm'))]
    refs=Counter(re.findall(r'R_X86_64_PLT32\s+('+PREFIX+r'\w+)-0x4',text))
    expected=doc(W/(stage+'-lower-result.json'))['core_calls']
    assert set(refs)==set(expected)
    assert all(0<refs[k]<=expected[k] for k in expected)
    if stage=='transport':
        exact=dict(expected);exact[PREFIX+'dense_divide8']-=1
        assert dict(refs)==exact
    symbols=subprocess.check_output(['/usr/bin/nm','--defined-only',str(path)],text=True)
    assert 'sirius_register_'+stage+'_packet' in symbols
    return dict(phase=stage+'-native',output=dict(path=str(path.relative_to(R)),**pin(path)),
        instructions=len(instructions),emitted_fused_instructions=[],core_relocations=dict(refs),
        IR_core_calls=expected,defined_symbols=symbols,
        native_call_site_decreases={k:expected[k]-refs[k] for k in expected if expected[k]!=refs[k]},
        scope='Compiler/object evidence only. IR call multiplicities remain exact by full-module normalized roundtrip. Any shared machine-call tail needs explicit CFG review before adoption; fewer static relocations alone never means fewer logical evaluations.')


phase,=sys.argv[1:]
if phase=='prepare':
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    failed=doc(W/'transport-native/owner.json')
    assert failed['status']=='terminal' and not failed['passed'] and failed['child_exit']==1
    assert failed['source_unchanged'] and failed['whole_inputs_unchanged'] and not failed['remaining_owned_processes']
    assert pin(O/'transport.o')==dict(bytes=3158032,sha256='00bd80cb8d795a78e8ed91087fe4b61dfd4cf94b3558c6c426cb605433ffbed8')
    plan=doc(W/'plan.json');inputs=dict(plan['inputs'])
    for p in [Path(__file__).resolve(),W/'transport-native/owner.json',W/'transport-native/stderr.log',O/'transport.o',W/'transport-lower-result.json',O/'transport.ll']:
        inputs[str(p.relative_to(R))]=pin(p)
    write(W/'resume-plan.json',dict(inputs=inputs,original_receipt_preserved=True,transport_recompiled=False))
    actions=doc(W/'actions.json')
    for name in ['audit-transport-existing','endpoint-lower-resume','endpoint-native-resume']:
        assert name not in actions
        extra=[]
        if name=='endpoint-lower-resume':extra=[O/'register_value_abi.so']
        if name=='endpoint-native-resume':extra=[O/'endpoint.bc',W/'endpoint-lower-result.json',W/'endpoint-lower-resume/owner.json']
        actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
            limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=['/usr/bin/python3','-B',str(Path(__file__).resolve()),name],
            frozen_inputs=sorted(inputs)+[str((W/'resume-plan.json').relative_to(R))]+[str(p.relative_to(R)) for p in extra],
            scope='Reuse completed Transport codegen after inspecting a static-call-count receipt failure. Compile only remaining Endpoint once. No changed FP/arithmetic/workloads/work limits or numerical/frame qualification.')
    (W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
    print(json.dumps(dict(prepared=['audit-transport-existing','endpoint-lower-resume','endpoint-native-resume'],transport_recompiled=False)))
else:
    for name,expected in doc(W/'resume-plan.json')['inputs'].items():assert pin(R/name)==expected,name
    if phase=='audit-transport-existing':
        result=native('transport',O/'transport.o')
        result['original_failed_owner']=pin(W/'transport-native/owner.json')
        write(W/'transport-native-result.json',result)
    elif phase=='endpoint-lower-resume':
        subprocess.run([sys.executable,'-B',str(W/'build.py'),'endpoint-lower'],check=True)
    elif phase=='endpoint-native-resume':
        owner=doc(W/'endpoint-lower-resume/owner.json');assert owner['passed'] and owner['status']=='terminal'
        lower=doc(W/'endpoint-lower-result.json')['output'];assert pin(R/lower['path'])=={k:lower[k] for k in ('bytes','sha256')}
        output=O/'endpoint.o';assert not output.exists()
        subprocess.run(doc(W/'plan.json')['commands']['endpoint-native'],check=True)
        result=native('endpoint',output)
        write(W/'endpoint-native-result.json',result)
    else:raise ValueError(phase)
    print(json.dumps(dict(phase=phase,compiler_repeated=False)))
