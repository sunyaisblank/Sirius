"""Build one fixed private-primitive comparison from exact source prefixes."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R=Path.cwd()
W=Path(__file__).resolve().parent
B=R/'bin/linux-gcc/native-cpu-exact-product-zero-residual'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())


def write(path,value):
    with path.open('x') as stream:stream.write(json.dumps(value,indent=2)+'\n')


def terminal(phase):
    value=doc(W/phase/'owner.json')
    assert value['status']=='terminal' and value['passed'] and value['child_exit']==0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    return value


phase,=sys.argv[1:]
proof=doc(W/'derivation.json')
for name,identity in proof['source_inputs'].items():assert pin(R/name)==identity,name
for name,identity in proof['derived'].items():assert pin(W/name)==identity,name
if phase in ('prepare-build','prepare-run'):
    assert not subprocess.check_output(['git','status','--porcelain']).strip()
    actions=doc(W/'actions.json');inputs=set(actions['build']['frozen_inputs'])
    inputs|={str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    if phase=='prepare-build':
        terminal('build')
        assert pin(B/'optional.o')==doc(W/'compile-result.json')['optional_object']
        headers={}
        for label,path in (('original',R/'out/native-cpu-value-fact-stores/retained8.cpp'),('candidate',W/'retained8.cpp')):
            source=path.read_text();start=source.index('namespace {')
            anchor='__attribute__((always_inline)) inline Triple finish('
            assert source.count(anchor)==1
            body=source[start:source.index(anchor)]
            assert body.count('namespace {')==1
            header=body.replace('namespace {','namespace '+label+' {',1)+'}\n'
            target=B/(label+'-prefix.hpp')
            with target.open('x') as stream:stream.write(header)
            assert header.replace('namespace '+label+' {','namespace {',1)[:-2]==body
            headers[str(target.relative_to(R))]=pin(target)
        args=list(proof['compile_argv']);args=args[:args.index('-c')]
        primitive=args+['-I',str(B),'-c',str(W/'primitive.cpp'),'-o',str(B/'primitive.o')]
        baseline=list(args);baseline[baseline.index('-march=x86-64-v4')]='-march=x86-64';baseline.remove('-mprefer-vector-width=256')
        driver=baseline+['-c',str(W/'driver.cpp'),'-o',str(B/'primitive-driver.o')]
        link=[args[0],str(B/'primitive-driver.o'),str(B/'primitive.o'),'-Wl,--no-relax','-o',str(B/'primitive-check')]
        plan=dict(headers=headers,compile_primitive=primitive,compile_driver=driver,link=link,
            original_source=pin(R/'out/native-cpu-value-fact-stores/retained8.cpp'),candidate_source=pin(W/'retained8.cpp'),
            exact_prefixes_through_product32=True,rest_of_arithmetic_not_in_fixture=True)
        write(W/'primitive-plan.json',plan);inputs|=set(headers)|{str((W/'primitive-plan.json').relative_to(R))}
        name,argv='primitive-build',['/usr/bin/python3','-B',str(W/'primitive.py'),'build']
    else:
        terminal('primitive-build')
        inputs|={str((B/n).relative_to(R)) for n in ('primitive.o','primitive-driver.o','primitive-check')}
        inputs|=set(doc(W/'primitive-plan.json')['headers'])
        name,argv='primitive-run',[str(B/'primitive-check')]
    assert name not in actions and not (W/name).exists()
    actions[name]=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        kind='diagnostic',expected_exit=0,software_provider=False,
        limits=dict(seconds=300,rss_kib=4194304,reserve_kib=2097152),argv=argv,frozen_inputs=sorted(inputs),
        scope='One fixed24-packet private-product comparison with exact original source prefix,19 independent integer-dyadic centers and signed-zero expectations, strict controls/sticky-state/canary/input checks. No random corpus, benchmark, provider/frame or scientific/release acceptance.')
    (W/'actions.json').write_text(json.dumps(actions,indent=2)+'\n')
    print(json.dumps(dict(prepared=name,executed=False)))
elif phase=='build':
    plan=doc(W/'primitive-plan.json')
    for name,identity in plan['headers'].items():assert pin(R/name)==identity,name
    for key in ('compile_primitive','compile_driver','link'):subprocess.run(plan[key],check=True)
    write(W/'primitive-build-result.json',dict(passed=True,plan=pin(W/'primitive-plan.json'),
        products={n:pin(B/n) for n in ('primitive.o','primitive-driver.o','primitive-check')}))
elif phase=='read':
    terminal('primitive-run')
    result=json.loads((W/'primitive-run/stdout.log').read_text())
    for key,expected in dict(passed=True,case_count=24,candidate_calls=24,historical_primitive_calls=24,
            compared_hi_low_words=384,compared_canary_words=384,independent_center_cases=19,denormal_reference_witness=True).items():
        assert result[key]==expected,(key,result[key])
    assert len(result['cases'])==24 and all(row['passed'] for row in result['cases'])
    assert sum(row['independent'] for row in result['cases'])==19
    write(W/'primitive-compatibility.json',dict(passed=True,result=result,
        products={n:pin(B/n) for n in ('primitive.o','primitive-driver.o','primitive-check')},
        limitations='Fixed private primitive source comparison only, including explicit signed zeros and all MXCSR sticky bits. Actual optional object code/retention dependency must be checked separately; no fast-branch PC/hit-rate, performance/provider/frame/scientific/release claim.'))
    print(json.dumps(dict(passed=True,cases=24,independent_centers=19,exact_zero_signs_and_MXCSR=True)))
else:
    raise ValueError(phase)
