"""Reuse optimized AU IR; change its private ABI, then compile each object once."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import shlex
import subprocess
import sys

R=Path.cwd()
W=Path(__file__).resolve().parent
V=R/'out/native-cpu-strict-ir-linking'
B=R/'bin/linux-gcc/native-cpu-strict-ir-linking'
D=B/'toolchain/bin'
O=R/'bin/linux-gcc/native-cpu-register-value'
INPUTS=dict(math='shared-math.O2.bc',dp='shared-dp.O2.bc',transport='shared-transport.O2.bc',endpoint='shared-endpoint.cleanup.O2.bc')
SUFFIXES=['dense_add8','dense_multiply8','dense_divide8','dense_sqrt8','tagged_add8','tagged_subtract8','tagged_multiply8','tagged_divide8','tagged_sqrt8']
PREFIX='sirius_register_value_'
COUNTS=dict(math=0,dp=1581,transport=2701,endpoint=3959)


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())


def write(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2);stream.write('\n')


phase,=sys.argv[1:]
if phase=='prepare':
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    assert not O.exists()
    O.mkdir()
    config=shlex.split(subprocess.check_output([str(D/'llvm-config'),'--cxxflags'],text=True))
    commands=dict(plugin=[str(D/'clang++'),*config,'-O2','-fPIC','-shared',str(V/'register_value_abi.cpp'),str(B/'sha256-only.a'),'-Wl,--exclude-libs,ALL','-o',str(O/'register_value_abi.so')])
    outputs=dict(plugin='register_value_abi.so')
    for stage,filename in INPUTS.items():
        lower,native=stage+'-lower',stage+'-native'
        outputs.update({lower:stage+'.bc',native:stage+'.o'})
        commands[lower]=[str(D/'opt'),'-load-pass-plugin='+str(O/'register_value_abi.so'),'-passes=sirius-register-'+stage+',verify','-verify-each',str(B/filename),'-o',str(O/outputs[lower])]
        commands[native]=[str(D/'llc'),'-O2','-filetype=obj','-relocation-model=pic','-mcpu='+('x86-64-v4' if stage=='math' else 'skylake-avx512')]
        if stage=='math':
            commands[native]+=['-print-after=finalize-isel','-filter-print-funcs='+','.join(PREFIX+s for s in SUFFIXES)]
        commands[native]+=[str(O/outputs[lower]),'-o',str(O/outputs[native])]
    paths=[*W.glob('*.py'),V/'register_value_abi.cpp',V/'toolchain.json',B/'sha256-only.a',
           R/'attestations/qualification/844dbeda/native-cpu-shared-value-original-frame-attempt/MANIFEST.json']
    paths+=[B/name for name in INPUTS.values()]
    for name in ('opt','llc','clang++','llvm-config','llvm-dis'):
        row=json.loads((V/'toolchain.json').read_text())['tools'][name]
        p=R/row['resolved'];assert pin(p)=={k:row[k] for k in ('bytes','sha256')};paths.append(p)
    inputs={str(p.relative_to(R)):pin(p) for p in paths}
    write(W/'plan.json',dict(commands=commands,outputs=outputs,inputs=inputs,
        convention='LLVM23 X86_RegCall(92), fixed across all4 separately compiled modules',
        scope='ABI/name setters only with exact normalized full-module roundtrip. No LLVM IR optimization, altered FP effects/arithmetic, workload, renderer or acceptance claim.'))
    actions={}
    for name in commands:
        deps=[] if name=='plugin' else ['register_value_abi.so']
        if name.endswith('-native'):deps.append(outputs[name.replace('-native','-lower')])
        actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
            limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),
            argv=['/usr/bin/python3','-B',str(Path(__file__).resolve()),name],
            frozen_inputs=sorted(inputs)+[str((W/'plan.json').relative_to(R))]+[str((O/name).relative_to(R)) for name in deps],
            scope='Compile one changed private ABI using original optimized AU modules. All strict FP/source instructions preserved exactly; original target feature closure and actual no-FMA gate retained. No numerical, benefit, frame or release acceptance.')
    write(W/'actions.json',actions)
    print(json.dumps(dict(prepared=list(commands),optimization_pipeline_repeated=False,workload_run=False)))
else:
    plan=json.loads((W/'plan.json').read_text());assert phase in plan['commands']
    for name,expected in plan['inputs'].items():assert pin(R/name)==expected,name
    predecessor=None if phase=='plugin' else ('plugin' if phase.endswith('-lower') else phase.replace('-native','-lower'))
    if predecessor:
        owner=json.loads((W/predecessor/'owner.json').read_text());assert owner['passed'] and owner['status']=='terminal'
        product=json.loads((W/(predecessor+'-result.json')).read_text())['output']
        assert pin(R/product['path'])=={k:product[k] for k in ('bytes','sha256')}
    output=O/plan['outputs'][phase];assert not output.exists()
    subprocess.run(plan['commands'][phase],check=True)
    result=dict(phase=phase,output=dict(path=str(output.relative_to(R)),**pin(output)),argv=plan['commands'][phase])
    if phase.endswith('-lower'):
        stage=phase.removesuffix('-lower');ll=output.with_suffix('.ll')
        subprocess.run([str(D/'llvm-dis'),str(output),'-o',str(ll)],check=True)
        source=ll.read_text()
        assert 'sirius_shared_value_' not in source and 'sirius_scope_'+stage+'_packet' not in source
        definitions=re.findall(r'^define .*$',source,re.M)
        calls=Counter(re.findall(r'\bcall\s+.*?@('+PREFIX+r'[^ (]+)\(',source))
        assert sum(calls.values())==COUNTS[stage]
        interfaces=re.findall(r'^(?:define|declare) .*@'+PREFIX+r'.*$',source,re.M)
        assert len(interfaces)==(8 if stage=='dp' else 9)
        assert all('x86_regcallcc ' in row for row in interfaces)
        core_calls=[row for row in source.splitlines() if re.search(r'\bcall\s+.*?@'+PREFIX,row)]
        assert all('x86_regcallcc ' in row for row in core_calls)
        attrs=dict(re.findall(r'^attributes #(\d+) = (.*)$',source,re.M))
        assert all('strictfp' in attrs[re.search(r'#(\d+)',row)[1]] for row in definitions)
        fast=re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',source,re.M)
        assert not fast
        result.update(ir=dict(path=str(ll.relative_to(R)),**pin(ll)),core_calls=dict(calls),interfaces=interfaces,definitions=len(definitions),normalized_roundtrip='exact; pass aborts on any non-ABI difference')
    elif phase.endswith('-native'):
        disassembly=subprocess.check_output(['/usr/bin/objdump','-drw',str(output)],text=True)
        instructions=re.findall(r'^\s*[0-9a-f]+:\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9.]*)',disassembly,re.M)
        fused=[x for x in instructions if x.startswith(('vfm','vfnm'))]
        assert instructions and not fused,fused
        refs=Counter(re.findall(r'R_X86_64_PLT32\s+('+PREFIX+r'\w+)-0x4',disassembly))
        stage=phase.removesuffix('-native');assert sum(refs.values())==COUNTS[stage]
        symbols=subprocess.check_output(['/usr/bin/nm','--defined-only',str(output)],text=True)
        if stage=='math':assert sum(PREFIX+s in symbols for s in SUFFIXES)==9
        else:assert 'sirius_register_'+stage+'_packet' in symbols
        result.update(instructions=len(instructions),emitted_fused_instructions=fused,core_relocations=dict(refs),defined_symbols=symbols,
            scope='Static object/relocation evidence only. Register assignment must be checked against math-native/stderr.log MIR; no speed claim from naming convention.')
    write(W/(phase+'-result.json'),result)
    print(json.dumps({k:v for k,v in result.items() if k in ('phase','output','definitions','instructions','core_calls','core_relocations')}))
