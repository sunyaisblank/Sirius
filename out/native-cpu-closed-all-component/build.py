"""One baseline adapter compile and11-part relocation link; no workload."""
from pathlib import Path
import json
import subprocess
import sys
from component import ingredients,doc,pin,checked,CORES,RAW_NAMES

R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/native-cpu-closed-all-component'

def write(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2);stream.write('\n')

phase,=sys.argv[1:]
if phase=='prepare':
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    value=ingredients();assert not B.exists() and not (W/'actions.json').exists()
    compile_argv=value['derivation']['compile_argv']
    merge_argv=['/usr/bin/ld','-r','-o',str(B/'rows_adapter.o'),str(B/'adapter-only.o')]+[str(R/row['path']) for row in value['products']]
    inputs=dict(value['inputs'])
    for p in [*W.glob('*.py'),W/'rows_adapter.cpp',W/'derivation.json',Path(compile_argv[0]).resolve(),
              Path('/usr/bin/ld').resolve(),Path('/usr/bin/nm').resolve(),Path('/usr/bin/readelf').resolve()]:
        inputs[str(p.relative_to(R)) if p.is_relative_to(R) else str(p)]=pin(p)
    write(W/'build-plan.json',dict(inputs=inputs,products=value['products'],compile_argv=compile_argv,merge_argv=merge_argv,
        new_native_compilations=1,native_calls=value['native_calls'],numerical_qualification_transferred=False,
        scope='Compile only baseline adapter after exact two private target-name substitutions; ld-r11 constituents. Old public source prefix/guards/FPscope/header/scalar/ISA fallbacks unchanged; native T/E counts bind actual backend evidence. Static interface only.'))
    write(W/'actions.json',{'build':dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        kind='diagnostic',expected_exit=0,software_provider=False,
        limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),
        argv=['/usr/bin/python3','-B',str(W/'build.py'),'build'],
        frozen_inputs=sorted(inputs)+['out/native-cpu-closed-all-component/build-plan.json'],
        scope='One baseline adapter compile plus11-part link/static symbols; no raw/math/provider/workload compilation or numerical acceptance transfer. Fresh actual original C8 remains required.')})
elif phase=='build':
    value=ingredients();plan=doc(W/'build-plan.json')
    assert plan['products']==value['products'] and plan['native_calls']==value['native_calls']
    for name,expected in plan['inputs'].items():checked(R/name,expected)
    assert not B.exists();B.mkdir()
    subprocess.run(plan['compile_argv'],check=True);subprocess.run(plan['merge_argv'],check=True)
    output=B/'rows_adapter.o'
    undefined=subprocess.check_output(['/usr/bin/nm','-u',str(output)],text=True)
    assert not any(prefix in undefined for prefix in ('sirius_shared_value_','sirius_local_value_','sirius_scope_','sirius_closed_'))
    listing=subprocess.check_output(['/usr/bin/readelf','-Ws',str(output)],text=True)
    cores,closed=[],[]
    for line in listing.splitlines():
        fields=line.split()
        if len(fields)!=8:continue
        name=fields[-1]
        assert not name.startswith(('sirius_shared_value_','sirius_scope_')),name
        if name in CORES or name in RAW_NAMES:
            size=int(fields[2],16 if fields[2].startswith('0x') else 10)
            assert fields[3:6]==['FUNC','GLOBAL','HIDDEN'] and fields[6]!='UND' and size>0
            (cores if name in CORES else closed).append(name)
    assert set(cores)==CORES and len(cores)==9 and set(closed)==set(RAW_NAMES) and len(closed)==3
    adapter=dict(path=str((B/'adapter-only.o').relative_to(R)),**pin(B/'adapter-only.o'))
    write(W/'build-result.json',dict(passed=True,compiled_interface_passed=True,
        output=dict(path=str(output.relative_to(R)),**pin(output)),adapter=adapter,parts=[adapter]+value['products'],
        core_definitions=sorted(cores),closed_raw_exports=list(RAW_NAMES),all_raw_scope_and_core_references_bound=True,
        new_native_compilations=1,numerical_qualification_transferred=False,remaining_external_symbols=undefined,
        scope='Actual baseline adapter compile/11-part static binding only; no scientific/numerical/runtime/benefit/frame acceptance transfer.'))
else:raise ValueError(phase)
