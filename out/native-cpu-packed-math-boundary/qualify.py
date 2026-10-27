"""One changed raw compilation and two existing complete-bank fixtures."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
import zipfile

R=Path.cwd();W=Path(__file__).resolve().parent
B=R/'bin/linux-gcc/native-cpu-packed-math-boundary';O=B/'objects'
sys.dont_write_bytecode=True
sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'))
sys.path.insert(0,str(R/'out/native-cpu-v4-raw-packets'))
from fixture_elf import associate


def doc(p):return json.loads(p.read_text())


def pin(p):
    with p.open('rb') as f:return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(f,'sha256').hexdigest())


def named(p):return dict(path=str(p.relative_to(R)),**pin(p))


def write(p,j):
    with p.open('x') as f:f.write(json.dumps(j,indent=2)+'\n')


def terminal(phase):
    o=doc(W/phase/'owner.json')
    assert o['status']=='terminal' and o['passed'] and o['source_unchanged'] and o['whole_inputs_unchanged']
    assert not o['remaining_owned_processes'] and not o['cleanup_errors']
    return o


phase,=sys.argv[1:];proof=doc(W/'derivation.json')
for key in ('source_inputs','derived'):
    for name,identity in proof[key].items():assert pin(R/name)==identity,name
math_path=R/'attestations/qualification/ca334e6e/native-cpu-exact-product-zero-residual-finite-compatibility/MANIFEST.json'
assert pin(math_path)['sha256']=='44d228597b02d724ad7a5b0ebfe8d69073bf306977bacd3a93013d1462213a66'
math=doc(math_path);assert math['finite_component_passed'] and not math['provider_adopted']
optional=R/'bin/linux-gcc/native-cpu-exact-product-zero-residual/optional.o'
assert pin(optional)==math['products'][str(optional.relative_to(R))]
sealed_path=R/'attestations/qualification/984ba457/native-cpu-v4-raw-packets-finite-compatibility/MANIFEST.json'
sealed=doc(sealed_path);archive=R/sealed['archive']['path']
assert pin(archive)=={k:sealed['archive'][k] for k in ('bytes','sha256')}

if phase.startswith('prepare-'):
    step=phase.removeprefix('prepare-');assert step in ('build','link','run')
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    inputs=set(proof['source_inputs'])|set(proof['derived'])|{str(math_path.relative_to(R)),str(optional.relative_to(R)),str(sealed_path.relative_to(R)),str(archive.relative_to(R))}
    inputs|={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    inputs|={'out/native-cpu-field-compatibility/elf_proof.py','out/native-cpu-accumulate-compatibility/elf_proof.py','out/native-cpu-v4-raw-packets/fixture_elf.py','out/native-cpu-row-ownership/cold-frame/run_owned.py','out/native-cpu-row-ownership/cold-frame/owned_guard.py','out/native-cpu-endpoint-header-provider/provider/source-contract.json'}
    for row in doc(R/'out/native-cpu-endpoint-header-provider/provider/source-contract.json')['native_objects'][:7]:inputs.add(row['path'])
    inputs|={'bin/linux-gcc/native-cpu-v4-raw-packets/objects/rows_adapter.o','bin/linux-gcc/native-cpu-v4-arithmetic/dispatch.o','bin/linux-gcc/native-cpu-v4-arithmetic/v3.o','bin/linux-gcc/native-cpu-accumulate-compatibility/baseline.o'}
    for name in ('/usr/bin/ld','/usr/bin/objdump','/usr/bin/nm','/usr/bin/ldd','/usr/bin/x86_64-linux-gnu-g++-14'):inputs.add(str(Path(name).resolve(strict=True)))
    actions={} if step=='build' else doc(W/'actions.json')
    if step!='build':
        terminal('build');assert doc(W/'code-proof.json')['passed']
        inputs|={str(p.relative_to(R)) for p in O.glob('*.o')}
    if step=='run':
        terminal('link');inputs|={str((B/(s+'-check')).relative_to(R)) for s in ('transport','endpoint')}
        inputs|=set(doc(W/'link-result.json')['runtime_libraries'])
    assert step not in actions and not (W/step).exists()
    actions[step]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),kind='diagnostic',expected_exit=0,software_provider=False,
        limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=['/usr/bin/python3','-B',str(W/'qualify.py'),step],frozen_inputs=sorted(inputs),
        scope='Changed private576B math-boundary raw3 with exact original tuple order/masks/FP and reviewed bounded zero-residual math. Source/build and two original finite fixtures only. No random cases, timing benchmark/provider/frame/scientific/release admission.')
    (W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n');print(json.dumps(dict(prepared=step,executed=False)))
elif phase=='build':
    O.mkdir(exist_ok=False)
    commands=proof['raw_compilations'];assert len(commands)==3
    for row in commands:
        assert not Path(row['output']).exists();subprocess.run(row['argv'],check=True)
    restored=[]
    with zipfile.ZipFile(archive) as z:
        for name in ('cpu_gate.o','transport-driver.o','endpoint-driver.o'):
            member='bin/linux-gcc/native-cpu-v4-raw-packets/objects/'+name;data=z.read(member)
            assert dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())==sealed['members'][member]
            destination=O/name
            with destination.open('xb') as f:f.write(data)
            restored.append(dict(path=str(destination.relative_to(R)),**pin(destination),archive_member=member))
    write(W/'build-result.json',dict(passed=True,commands=commands,new_raw_compilations=3,new_CXX_compilations=0,restored_fixture_objects=restored,objects=[named(O/(s+'.o')) for s in ('transport','endpoint','dopri_phase')]))
elif phase=='link':
    assert doc(W/'code-proof.json')['passed']
    old=doc(R/'out/native-cpu-endpoint-header-provider/provider/source-contract.json')
    objects=[R/x['path'] for x in old['native_objects'][:7]]
    objects += [R/'bin/linux-gcc/native-cpu-v4-arithmetic'/n for n in ('dispatch.o','v3.o')]+[optional]
    objects += [O/(s+'.o') for s in ('transport','endpoint','dopri_phase')]+[O/'cpu_gate.o']
    results=[];libraries={}
    for stage in ('transport','endpoint'):
        selected=[O/(stage+'-driver.o'),*objects]
        if stage=='transport':selected += [R/'bin/linux-gcc/native-cpu-accumulate-compatibility/baseline.o',R/'bin/linux-gcc/native-cpu-v4-raw-packets/objects/rows_adapter.o']
        exe=B/(stage+'-check');assert not exe.exists()
        argv=['/usr/bin/x86_64-linux-gnu-g++-14',*map(str,selected),'-Wl,--no-relax','-lm','-o',str(exe)]
        subprocess.run(argv,check=True)
        joined=associate(exe,[dict(path=str(p),**pin(p)) for p in selected]);write(W/(stage+'-link-proof.json'),joined)
        resolution=subprocess.check_output(['ldd','-r',str(exe)],text=True,stderr=subprocess.STDOUT)
        assert 'not found' not in resolution and 'undefined symbol' not in resolution
        for name in re.findall(r'(/[^\s()]+)',resolution):
            p=Path(name)
            if p.is_file():p=p.resolve(strict=True);libraries[str(p)]=pin(p)
        results.append(dict(stage=stage,argv=argv,executable=named(exe)))
    write(W/'link-result.json',dict(passed=True,executables=results,runtime_libraries=libraries,new_CXX_compilations=0))
elif phase=='run':
    results=[]
    for stage in ('transport','endpoint'):
        argv=[str(B/(stage+'-check'))]
        if stage=='transport':folder=W/'run/transport';folder.mkdir();argv.append(str(folder))
        result=subprocess.run(argv,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        comparison=json.loads(result.stdout);assert comparison['passed'] and not result.stderr
        write(W/(stage+'-output.json'),dict(stdout=result.stdout,stderr=result.stderr));results.append(dict(stage=stage,result=comparison))
        print(result.stdout,end='')
    write(W/'compatibility.json',dict(passed=True,fixtures=results,limitations='Original Transport38 and focused Endpoint complete-word fixtures only, coupled with source-qualified new arithmetic. DopriPhase compiled/linked but not executed; no performance, provider/C8/all8/ABBA/frame/scientific/release admission.'))
else:raise ValueError(phase)
