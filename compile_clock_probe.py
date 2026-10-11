"""Root-owned finite existing-toolchain compiler/validator admission, no GPU work."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,os,signal,subprocess,time
ROOT=Path.cwd();W=Path(__file__).resolve().parent
compiler=Path('/opt/slang/bin/slangc');library=Path('/opt/slang/lib/libslang-compiler.so.0.2026.12.0.1')
tools=ROOT/'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def birth(pid):
 try:return Path('/proc',str(pid),'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
assert seal(compiler)=={'bytes':60464,'sha256':'26d7fd442248bc702dd07c07757d628e5bf20c17a95359aa54191ab19da4ff7b'}
assert seal(library)=={'bytes':29840880,'sha256':'c32802f3627b7956ab5f5413eec69db8f1b7da9f3166e4e50dcbefbfcbef4b72'}
paths=[Path(__file__),W/'clock_probe.slang',compiler,library,tools/'spirv-val',tools/'spirv-dis',ROOT/'scripts/build-retained-kernels.py']
record={'input_seals':{str(p):seal(p) for p in paths},'seconds_per_command':60,'rss_bytes':536870912,'runs':[],'scope':'Compiler and SPIR-V validator admission only; no logical device, pipeline, dispatch or clock timing verdict.'}
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
  (W/'compiler-admission.json').write_text(json.dumps(record,indent=2)+'\n')
 assert r['accepted'],r
run('probe-slang',[compiler,W/'clock_probe.slang','-O0','-target','spirv','-profile','spirv_1_5','-entry','ComputeMain','-stage','compute','-denorm-mode-fp32','preserve','-o',W/'clock_probe.spv'])
run('probe-validate',[tools/'spirv-val','--target-env','vulkan1.2',W/'clock_probe.spv'])
run('probe-disassemble',[tools/'spirv-dis',W/'clock_probe.spv','-o',W/'clock_probe.spvasm'])
text=(W/'clock_probe.spvasm').read_text()
assert 'OpCapability ShaderClockKHR' in text and 'OpExtension "SPV_KHR_shader_clock"' in text and text.count('OpReadClockKHR')==1
assert 'OpTypeInt 64' not in text and 'OpCapability Int64' not in text
import re
u32=re.search(r'(\S+) = OpTypeInt 32 0',text)[1];uvec2=re.search(r'(\S+) = OpTypeVector '+re.escape(u32)+r' 2',text)[1]
scope=re.search(r'(\S+) = OpConstant '+re.escape(u32)+r' 1\b',text)[1]
assert re.search(r'OpReadClockKHR '+re.escape(uvec2)+' '+re.escape(scope)+r'\b',text)
assert {str(p):seal(p) for p in paths}==record['input_seals']
record.update(accepted=True,result={'scope':'Device(1)','type':'uint2','clock_reads':1,'shaderInt64_required':False},outputs={str(p):seal(p) for p in (W/'clock_probe.spv',W/'clock_probe.spvasm')},signals=signals)
(W/'compiler-admission.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'accepted':True,'result':record['result'],'runs':[{k:r[k] for k in ('name','elapsed_seconds','peak_rss_bytes','returncode')} for r in record['runs']]}))
