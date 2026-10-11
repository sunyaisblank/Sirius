"""Root-owned finite existing-toolchain compiler/validator admission, no GPU work."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,hashlib,json,os,re,signal,subprocess,time
ROOT=Path.cwd();W=Path(__file__).resolve().parent
compiler=Path('/opt/slang/bin/slangc');library=Path('/opt/slang/lib/libslang-compiler.so.0.2026.12.0.1')
tools=ROOT/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
assert seal(compiler)=={'bytes':60464,'sha256':'26d7fd442248bc702dd07c07757d628e5bf20c17a95359aa54191ab19da4ff7b'}
assert seal(library)=={'bytes':29840880,'sha256':'c32802f3627b7956ab5f5413eec69db8f1b7da9f3166e4e50dcbefbfcbef4b72'}
paths=[Path(__file__),W/'phase/retained_transport.slang',W/'phase/retained_endpoint.slang',W/'phase-preparation.json',W/'clock-features.json',compiler,library,tools/'spirv-val',tools/'spirv-dis',tools/'spirv-as',ROOT/'scripts/build-retained-kernels.py']
commands_file=ROOT/'attestations/review/b81e061/native-narrow-fma-premise/diagnostic/native-narrow-fma-premise/compiler-commands.json'
paths.append(commands_file)
paths += [W/'failed-phase-compiler-v0/phase-compiler.json',W/'phase/retained_transport.compiler.spv']
paths.extend((ROOT/'src/sirius/kernels').glob('retained_*.slang'))
record={'input_seals':{str(p):seal(p) for p in paths},'seconds_per_command':60,'rss_bytes':536870912,'runs':[],'scope':'Unadopted native-default complete-module compilation and SPIR-V inspection only; no logical device, pipeline, dispatch, driver-placement or clock timing verdict.'}
signals=[]
for s in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):signal.signal(s,lambda n,f:signals.append(n))
def run(name,args):
 r={'name':name,'arguments':list(map(str,args)),'errors':[],'outer_stop':False,'peak_rss_bytes':0};p=None;t=time.monotonic()
 try:
  with (W/(name+'.stdout')).open('xb') as out,(W/(name+'.stderr')).open('xb') as err:
   p=subprocess.Popen(r['arguments'],stdout=out,stderr=err,start_new_session=True);ticks=birth(p.pid);r.update(pid=p.pid,start_ticks=ticks)
   while p.poll() is None:
    try:rss=int(Path('/proc',str(p.pid),'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE');r['peak_rss_bytes']=max(r['peak_rss_bytes'],rss)
    except FileNotFoundError:pass
    if signals or time.monotonic()-t>60 or r['peak_rss_bytes']>536870912:raise RuntimeError('fixed signal/time/RSS guard triggered')
    try:p.wait(timeout=.05)
    except subprocess.TimeoutExpired:pass
   assert p.returncode==0
 except BaseException as e:r['errors'].append(repr(e));r['outer_stop']=True
 finally:
  if p:
   if p.poll() is None:
    assert birth(p.pid)==ticks and os.getpgid(p.pid)==p.pid;os.killpg(p.pid,signal.SIGKILL)
   r['returncode']=p.wait(timeout=10);r['owned_birth_absent']=birth(p.pid)!=ticks
   try:os.killpg(p.pid,0);r['owned_group_absent']=False
   except ProcessLookupError:r['owned_group_absent']=True
  r['elapsed_seconds']=time.monotonic()-t;r['accepted']=not r['errors'] and not signals and r.get('returncode')==0 and r.get('owned_birth_absent') and r.get('owned_group_absent');record['runs'].append(r)
  (W/'phase-compiler.json').write_text(json.dumps(record,indent=2)+'\n')
 assert r['accepted'],r

features=json.loads((W/'clock-features.json').read_text())
assert features['clock_features']['shaderDeviceClock'] and features['logical_devices_created']==0
# Execute only the current builder's pure native-controls function.
builder=ROOT/'scripts/build-retained-kernels.py';tree=ast.parse(builder.read_text())
node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='native_controls')
env={'re':re};exec(compile(ast.Module(body=[node],type_ignores=[]),str(builder),'exec'),env)
commands=json.loads(commands_file.read_text());products={}
for stem,expected_calls in [('retained_transport',4),('retained_endpoint',3)]:
 command=next(c['argv'] for c in commands if Path(c['argv'][0]).name=='slangc' and Path(c['argv'][1]).name==stem+'.slang' and 'baseline-'+stem+'.compiler.spv' in c['argv'][-1])
 args=list(command);args[1]=str(W/'phase'/(stem+'.slang'));args[args.index('-I')+1]=str(ROOT/'src/sirius/kernels');args[-1]=str(W/'phase'/(stem+'.compiler.spv'))
 raw=W/'phase'/(stem+'.compiler.spv');assembly=W/'phase'/(stem+'.spvasm');final=W/'phase'/(stem+'.spv')
 if stem=='retained_transport':
  prior=json.loads((W/'failed-phase-compiler-v0/phase-compiler.json').read_text())
  r=prior['runs'][0];assert r['accepted'] and r['arguments']==list(map(str,args))
  record['reused_transport_compiler_run']=r
 else:run(stem+'-slang',args)
 run(stem+'-corrected-disassemble',[tools/'spirv-dis',raw,'-o',assembly])
 text=assembly.read_text();text=env['native_controls'](text);assembly.write_text(text)
 run(stem+'-corrected-assemble',[tools/'spirv-as','--target-env','spv1.5',assembly,'-o',final])
 run(stem+'-corrected-validate',[tools/'spirv-val','--target-env','vulkan1.2',final])
 assert text.count('OpReadClockKHR')==1 and 'OpCapability ShaderClockKHR' in text
 assert 'OpTypeInt 64' not in text and 'OpCapability Int64' not in text
 assert 'OpExtension "SPV_KHR_shader_clock"' in text
 assert re.search(r'OpExecutionMode %\S+ LocalSize 64 1 1',text)
 name=re.search(r'OpName (%\S+) "ObservePhase"',text)[1]
 calls=len(re.findall(r'OpFunctionCall %\S+ '+re.escape(name)+r'\b',text))
 assert calls==expected_calls,(stem,calls,expected_calls)
 products[stem]={'spv':seal(final),'assembly':seal(assembly),'clock_instruction_count':1,'clock_helper_call_sites':calls,'dynamic_success_clock_slots':16 if stem=='retained_transport' else 4,'clock_capability':'ShaderClockKHR','clock_scope':'Device(1)','logical_binding':2,'local_size':[64,1,1]}
 assert re.search(r'OpDecorate %phaseObservations Binding 2',text)
assert {str(p):seal(p) for p in paths}==record['input_seals']
record.update(accepted=True,products=products,signals=signals)
(W/'phase-compiler.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'accepted':True,'products':products,'runs':[{k:r[k] for k in ('name','elapsed_seconds','peak_rss_bytes','returncode')} for r in record['runs']]}))
