from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os, signal, subprocess, sys, time, xml.etree.ElementTree as ET

r = Path.cwd(); w = Path(__file__).resolve().parent
def git(*args): return subprocess.check_output(['git', *args], text=True).strip()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p, x): p.write_text(json.dumps(x, indent=2)+'\n')
def birth(pid): return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
def members(group):
    result=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdecimal(): continue
        try: f=(p/'stat').read_text().rsplit(')',1)[1].split()
        except (FileNotFoundError,ProcessLookupError): continue
        if int(f[2])==group and int(f[3])==group: result.append({'pid':int(p.name),'start_ticks':f[19]})
    return result
def resources(group):
    rss=0
    for x in members(group):
        try: s=Path('/proc',str(x['pid']),'status').read_text()
        except (FileNotFoundError,ProcessLookupError): continue
        rss+=next((int(v.split()[1]) for v in s.splitlines() if v.startswith('VmRSS:')),0)
    return rss
def kill_group(group,action):
    try: os.killpg(group,action)
    except ProcessLookupError: pass

assert len(sys.argv)==3 and sys.argv[1] in ('candidate','baseline')
mode,label=sys.argv[1:]
head=git('rev-parse','HEAD')
assert not git('status','--porcelain')
assert not {k:v for k,v in os.environ.items() if k.startswith('GTEST_') and v}
base=w/'baseline-8d5b89d'
if mode=='candidate':
    build=json.loads((w/'build-and-arrays.json').read_text())
    bindings=json.loads((w/'linux-readonly-payloads.json').read_text())
    assert build['source_revision']==head and build['returncode']==0 and build['all_36_unaffected_arrays_exact']
    assert bindings['source_revision']==head and bindings['pass_']
    header=r/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
    assert sha(header)==build['whole_header_sha256']==bindings['whole_header_sha256']
    artifact=r/'bin/linux-gcc/tests/backend/sirius_backend_tests'
    expected=bindings['actual_consumers']['bin/linux-gcc/tests/backend/sirius_backend_tests']
    source=head
    controls=[('dense','RetainedComputeTest.DenseSegmentsPreserveSmallCovariantArrivalDerivatives')]
    controls.append(('zero-arrival','RetainedComputeTest.ZeroFractionArrivalsPreserveProgramAuthorityAndCompleteRefusal'))
    if label.startswith('benefit-'): controls=[]
    controls.append(('sampler','RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates'))
else:
    source='8d5b89d15d8591db25093370368e12891a55bf03'
    bindings=json.loads((base/'linux-readonly-payloads.json').read_text())
    assert bindings['source_revision']==source and bindings['pass_']
    artifact=base/'sirius_backend_tests'
    expected=bindings['actual_consumers']['bin/linux-gcc/tests/backend/sirius_backend_tests']
    header=r/'attestations/software-vulkan/8d5b89d/retained-preparation-count/executed-inputs/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
    assert sha(header)==bindings['whole_header_sha256']
    controls=[('sampler','RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates')]
assert artifact.stat().st_size==expected['bytes'] and sha(artifact)==expected['sha256']
provider=[Path('/usr/share/vulkan/icd.d/lvp_icd.json'),Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so'),Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1').resolve()]
selected={'VK_DRIVER_FILES':str(provider[0]),'SIRIUS_VULKAN_DEVICE':'','MESA_SHADER_CACHE_DISABLE':'','LP_NATIVE_VECTOR_WIDTH':'','LP_PERF':'','GALLIVM_PERF':''}
files=[Path(__file__).resolve(),artifact,header,*provider,*sorted((r/'src/sirius/kernels').glob('retained*')),
       r/'scripts/build-retained-kernels.py',r/'tests/backend/retained_compute_test.cpp',r/'tests/backend/retained_dopri_test.cpp',
       r/'src/sirius/backend/retained_integrator.cpp',r/'src/sirius/backend/retained_compute.cpp',r/'src/sirius/backend/retained_trace_executor.cpp',
       *sorted((r/'tests/support/retained_transport').glob('*.h'))]
files+=([w/'source.json',w/'build-and-arrays.json',w/'linux-readonly-payloads.json',r/'bin/linux-gcc/generated/sirius/alignment_receipt.json'] if mode=='candidate'
        else [base/'source-binding.json',base/'linux-readonly-payloads.json',base/'build-and-arrays.json'])
seals={str(p):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in files}
receipts=[]
for name,case in controls:
    d=w/'controls'/(mode+'-'+label+'-'+name); d.mkdir(parents=True,exist_ok=False)
    args=[str(artifact),'--gtest_filter='+case,'--gtest_output=xml:'+str(d/'tests.xml')]
    owner={'execution_source_revision':source,'execution_source_tree':git('rev-parse',source+'^{tree}'),
           'checkout_revision':head,'checkout_clean':True,'mode':mode,'label':label,'case':case,'arguments':args,
           'whole_input_seals':seals,'selected_environment':selected,'timeout_seconds':180,'maximum_group_rss_mib':4096,
           'sample_cadence_seconds':.25,'controller_pid':os.getpid(),'controller_start_ticks':birth(os.getpid()),
           'peak_sampled_rss_kib':0,'stop_reason':None,'owner_error':None,'started_utc':datetime.now(timezone.utc).isoformat(),
           'scope':'Existing full Dense and/or original full moving-sampler GoogleTest case, original settings/counts/references/assertions. Baseline is preserved whole8d ELF, candidate is exact current clean build. Whole-case elapsed includes initialization, assertions, runtime observers and teardown; no isolated arithmetic, exclusive device, frame, native or fullqualification claim.'}
    env=os.environ.copy();env.update(selected)
    births={};loaded=set();started=time.monotonic()
    with (d/'stdout.log').open('wb') as out,(d/'stderr.log').open('wb') as err:
        child=subprocess.Popen(args,cwd=r,env=env,stdout=out,stderr=err,start_new_session=True)
        owner.update(child_pid=child.pid,child_start_ticks=birth(child.pid));births[child.pid]=owner['child_start_ticks']
        dump(d/'owner.json',owner)
        try:
            while child.poll() is None:
                births.update((x['pid'],x['start_ticks']) for x in members(child.pid))
                owner['peak_sampled_rss_kib']=max(owner['peak_sampled_rss_kib'],resources(child.pid))
                try:
                    for line in Path('/proc',str(child.pid),'maps').read_text().splitlines():
                        p=line.split()[-1]
                        if p.startswith('/') and ('libvulkan_lvp.so' in p or 'libvulkan.so' in p): loaded.add(str(Path(p).resolve()))
                except (FileNotFoundError,ProcessLookupError): pass
                if time.monotonic()-started>=180:owner['stop_reason']='case_time_bound'
                elif owner['peak_sampled_rss_kib']>4096*1024:owner['stop_reason']='case_rss_bound'
                dump(d/'owner.json',owner)
                if owner['stop_reason']:break
                time.sleep(.25)
        except BaseException as error:
            owner['owner_error']=repr(error);owner['stop_reason']=owner['stop_reason'] or 'observer_exception'
        finally:
            if child.poll() is None or members(child.pid):
                owner['stop_reason']=owner['stop_reason'] or 'remaining_owned_group'
                kill_group(child.pid,signal.SIGTERM)
                try:child.wait(timeout=10)
                except subprocess.TimeoutExpired:kill_group(child.pid,signal.SIGKILL);child.wait(timeout=10)
            child.wait()
    owner.update(returncode=child.returncode,elapsed_seconds=time.monotonic()-started,
                 observed_loaded_provider_modules=sorted(loaded),known_births=[{'pid':p,'start_ticks':s} for p,s in sorted(births.items())],
                 owned_group_remaining=members(child.pid),unchanged_inputs=all(p.stat().st_size==seals[str(p)]['bytes'] and sha(p)==seals[str(p)]['sha256'] for p in files),
                 ended_utc=datetime.now(timezone.utc).isoformat())
    xml=d/'tests.xml'
    owner['xml']=None
    if xml.exists():
        t=ET.parse(xml);root=t.getroot();cases=list(root.iter('testcase'))
        owner['xml']={'sha256':sha(xml),'tests':int(root.attrib.get('tests','0')),'failures':int(root.attrib.get('failures','0')),
                      'errors':int(root.attrib.get('errors','0')),'skipped':int(root.attrib.get('disabled','0'))+int(root.attrib.get('skipped','0')),
                      'cases':[{'name':x.attrib.get('classname','')+'.'+x.attrib.get('name',''),'attributes':x.attrib,'failure_elements':len(list(x.iter('failure'))),'skip_elements':len(list(x.iter('skipped')))} for x in cases]}
    x=owner['xml']
    owner['pass_']=child.returncode==0 and owner['stop_reason'] is None and owner['owner_error'] is None and not owner['owned_group_remaining'] and owner['unchanged_inputs'] and str(provider[1].resolve()) in loaded and x is not None and x['tests']==1 and x['failures']==x['errors']==x['skipped']==0 and len(x['cases'])==1 and x['cases'][0]['name']==case and x['cases'][0]['failure_elements']==x['cases'][0]['skip_elements']==0
    assert git('rev-parse','HEAD')==head and not git('status','--porcelain')
    dump(d/'owner.json',owner);receipts.append(owner)
    dump(w/(mode+'-'+label+'-sequence.json'),{'receipts':receipts,'pass_':all(q['pass_'] for q in receipts)})
    print(json.dumps({k:owner[k] for k in ('mode','label','case','returncode','elapsed_seconds','peak_sampled_rss_kib','stop_reason','pass_')}),flush=True)
    if not owner['pass_']:sys.exit(1)
