"""Bind one zero-dispatch ordinary driver query to the private peeling module."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
OLD=ROOT/'attestations/native-vulkan/451e2cb/native-division-phase-observer/files/placement'
ORDINARY=ROOT/'attestations/native-vulkan/3045067/ordinary-pipeline-generated-code'
BASELINE='2053ef79c90752d37b0075b8137cfb9b571b6820'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def git(*args):return subprocess.check_output(['git',*args],text=True,timeout=30).strip()
def key(p):return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p.resolve())
def one(s,old,new):
    assert s.count(old)==1,(old,s.count(old))
    return s.replace(old,new)
assert ROOT.resolve()==Path('/home/astra/.project/Sirius') and WORK==ROOT/'out/native-division-first-correction-screen'
assert git('rev-parse','HEAD')==BASELINE==git('rev-parse','@{u}') and not git('status','--porcelain')
compiler=read(WORK/'compiler.json');compiled=read(WORK/'compiled-review.json')
assert compiler['accepted'] and compiled['accepted'] and compiler['baseline_revision']==compiled['baseline_revision']==BASELINE
required={WORK/'compiler.json',WORK/'products/transport.spv',WORK/'products/transport.spvasm',WORK/'products/transport.compiler.spv',
          WORK/'products/original-transport.spvasm',WORK/'preparation.json',WORK/'source-review.json',WORK/'compiler-controller-review.json'}
required.update(Path(n) for n in compiler['input_seals'])
assert all(str(p) in compiled['reviewed_seals'] for p in required)
for path,value in compiled['reviewed_seals'].items():assert seal(Path(path))==value,path
assert seal(WORK/'products/transport.spv')==compiler['products']['spv']
prior=read(OLD/'source-review.json');assert prior['accepted']
target=WORK/'placement';assert not target.exists();target.mkdir()
records={}
for name in ('native_clock_placement.cs','placement-owner.ps1','placement-child.ps1','run_placement_query.py',
             'emergency_native_cleanup.ps1','run_native_retirement.py','verify_native_retirement.ps1'):
    source=OLD/name;assert seal(source)==prior['controller_seals'][name]
    text=source.read_text();changes=[]
    def replace(old,new):
        global text
        text=one(text,old,new);changes.append({'old':old,'new':new})
    if name=='native_clock_placement.cs':
        replace('// Unadopted native-default phase-observer ordinary driver-code query.\n// Explicit Device-clock feature and binding2; no dispatch or qualification claim.',
                '// Private native-default first-correction peeling ordinary driver-code query.\n// Original production feature chain and bindings0/1; no dispatch or qualification claim.')
        replace('    [StructLayout(LayoutKind.Sequential)] struct ClockFeatures {\n        public uint sType; public IntPtr pNext; public uint subgroupClock, deviceClock;\n    }\n','')
        replace('        Abi<ClockFeatures>(24);\n        if (Marshal.OffsetOf(typeof(ClockFeatures),"subgroupClock").ToInt32()!=16 || Marshal.OffsetOf(typeof(ClockFeatures),"deviceClock").ToInt32()!=20) throw new Exception("Clock feature ABI mismatch");\n','')
        replace(' || !advertised.Contains("VK_KHR_shader_clock")','')
        replace('            IntPtr shaderClock = Allocate(24, owned); Marshal.WriteInt32(shaderClock,1000181000);\n','')
        replace('IntPtr fma = Allocate(32, owned); Marshal.WriteInt32(fma, 1000579000); Marshal.WriteIntPtr(fma,8,shaderClock);',
                'IntPtr fma = Allocate(32, owned); Marshal.WriteInt32(fma, 1000579000);')
        replace('            uint subgroupClockSupported=U32(shaderClock,16),deviceClockSupported=U32(shaderClock,20);\n            if(subgroupClockSupported>1 || deviceClockSupported!=1 || U32(shaderClock,0)!=1000181000 || Marshal.ReadIntPtr(shaderClock,8)!=IntPtr.Zero) throw new Exception("Actual Device-clock feature unavailable or malformed");\n            Marshal.WriteInt32(shaderClock,16,0); Marshal.WriteInt32(shaderClock,20,1);\n','')
        replace('"VK_KHR_shader_fma", "VK_AMD_shader_info", "VK_KHR_shader_clock"','"VK_KHR_shader_fma", "VK_AMD_shader_info"')
        replace('IntPtr enabledExtensions = Allocate(24, owned);','IntPtr enabledExtensions = Allocate(16, owned);')
        replace('features = enabled, extensionCount = 3, extensions = enabledExtensions','features = enabled, extensionCount = 2, extensions = enabledExtensions')
        replace('var setInfo = new SetLayoutInfo { sType = 32, count = 3,','var setInfo = new SetLayoutInfo { sType = 32, count = 2,')
        replace(', new Binding { binding = 2, type = 7, count = 1, stages = 0x20 }','')
        replace('f1306ba2448bf0ce09bcbd55846fb8d6bbcddbc467d453a37e46fc3210f2cff4',seal(WORK/'products/transport.spv')['sha256'])
        replace('1c4687976bb4efb00d2c33c55a9690d0128591ab22188903894e47217959576d',seal(ORDINARY/'query-payloads/endpoint.spv')['sha256'])
        replace('"descriptor_bindings",new uint[]{0,1,2},"shader_device_clock_supported",deviceClockSupported,"shader_subgroup_clock_supported",subgroupClockSupported,"shader_device_clock_enabled",true,"shader_subgroup_clock_enabled",false',
                '"descriptor_bindings",new uint[]{0,1},"shader_device_clock_enabled",false,"shader_subgroup_clock_enabled",false')
        replace('One exact unadopted native-default division primitive observer module in a fresh process; Device-clock feature and descriptorbinding2 explicitly enabled with AMD tooling extension. Empty pipeline cache/flags0 ordinary compilation; driver coldness unverified. Static generated-code query only; no dispatch, clock execution/rate, timing benefit, numerical/scientific/release verdict.',
                'One exact private native-default first-correction peeling module in a fresh process; original production Float64/FMA32 feature chain and descriptorbindings0/1 with AMD tooling extension. Empty pipeline cache/flags0 ordinary compilation; driver coldness unverified. Static generated-code query only; no dispatch, execution timing benefit, numerical/scientific/release verdict.')
        assert 'shaderClock' not in text and 'VK_KHR_shader_clock' not in text
        assert 'for (int i = 0; i < 1; ++i)' in text
    elif name=='placement-owner.ps1':
        replace('Ordinary generated-code query for one fixed Transport diagnostic module in a fresh process with Device-clock enabled and binding2;',
                'Ordinary generated-code query for one fixed private first-correction peeling Transport module in a fresh process with original feature chain and bindings0/1;')
    elif name=='run_placement_query.py':
        replace("ROOT/'out/native-division-phase-observer/placement'","ROOT/'out/native-division-first-correction-screen/placement'")
        replace("WORK.parent/'prepare_division_placement.py'","WORK.parent/'prepare_query.py'")
    (target/name).write_text(text)
    records[name]={'basis':str(source.relative_to(ROOT)),'basis_seal':seal(source),'current_seal':seal(target/name),
                   'byte_identical':source.read_bytes()==(target/name).read_bytes(),'changes':changes}
for name,path in [('retained_transport.spv',WORK/'products/transport.spv'),('retained_endpoint.spv',ORDINARY/'query-payloads/endpoint.spv')]:
    (target/name).write_bytes(path.read_bytes());assert seal(target/name)==seal(path)
providers=read(OLD/'plan.json')['providers']
for provider in providers:assert seal(Path(provider['linux_path']))=={k:provider[k] for k in ('bytes','sha256')}
paths=set(Path(n) for n in compiled['reviewed_seals'])
paths.update((Path(__file__),WORK/'compiled-review.json',WORK/'compiler.json',WORK/'preparation.json',WORK/'source-review.json',
              OLD/'source-review.json',OLD/'plan.json',ORDINARY/'query-payloads/endpoint.spv'))
paths.update(p for p in target.iterdir() if p.is_file())
plan={'source_revision':BASELINE,'selected_stage':'Transport','input_seals':{key(p):seal(p) for p in sorted(paths)},'providers':providers,
      'purpose':'One exact private native-default first-correction peeling Transport module; ordinary static shader-info query, zero dispatch. Originalproduction Float64/FMA32, bindings0/1, emptycache/flags0 plus AMDqueryextension. Existing50s/512MiB child/45s intrinsicwatchdog/120sbridge/8MiBpayload guards unchanged. Actual compiler-specialization/footprint screen; no cost, runtime equivalence/benefit or qualification claim.'}
(target/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
(target/'preparation.json').write_text(json.dumps({'sources':records,'plan_seal':seal(target/'plan.json'),
      'scope':'Exact source preparation and original numerical feature/layout restoration only; no helper/native/device invocation.'},indent=2)+'\n')
assert git('rev-parse','HEAD')==BASELINE==git('rev-parse','@{u}') and not git('status','--porcelain')
print(json.dumps({'sources':{n:r['current_seal'] for n,r in records.items()},'plan':seal(target/'plan.json')}))
