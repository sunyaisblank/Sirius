"""Adapt the protected ordinary AMD query and reviewed owned bridge."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,hashlib,json
R=Path.cwd();W=Path(__file__).resolve().parent
head=json.loads((W/'plan.json').read_text())['source_revision']
def one(s,a,b):assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)
def save(n,s):p=W/n;assert not p.exists();p.write_text(s)
def seal(p):return {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
s=(W/'prior-native_endpoint_query.cs').read_text();s=one(s,'public static class NativeTransportFmaProbe {','public static class NativeClockPlacement {')
s=s.replace('// Baseline-only ordinary native pipeline generated-code observer.\n// Adds only the advertised AMD query extension; no dispatch or qualification claim.','// Unadopted native-default phase-observer ordinary driver-code query.\n// Explicit Device-clock feature and binding2; no dispatch or qualification claim.')
s=one(s,'    [StructLayout(LayoutKind.Sequential)] struct ApplicationInfo {','    [StructLayout(LayoutKind.Sequential)] struct ClockFeatures {\n        public uint sType; public IntPtr pNext; public uint subgroupClock, deviceClock;\n    }\n    [StructLayout(LayoutKind.Sequential)] struct ApplicationInfo {')
s=one(s,'        Abi<ApplicationInfo>(48);','        Abi<ClockFeatures>(24);\n        if (Marshal.OffsetOf(typeof(ClockFeatures),"subgroupClock").ToInt32()!=16 || Marshal.OffsetOf(typeof(ClockFeatures),"deviceClock").ToInt32()!=20) throw new Exception("Clock feature ABI mismatch");\n        Abi<ApplicationInfo>(48);')
anchor='    public static object Observe(string loaderPath, string transportPath, string endpointPath, string outputDirectory) {'
wrapper='''    public static object Observe(string loaderPath, string transportPath, string endpointPath, string outputDirectory, int parentPid, string parentBirth) {
        bool done=false;var timer=System.Diagnostics.Stopwatch.StartNew();
        var watchdog=new System.Threading.Thread(delegate(){
            while(!System.Threading.Volatile.Read(ref done)){
                try {using(var parent=System.Diagnostics.Process.GetProcessById(parentPid)){
                    if(parent.StartTime.ToUniversalTime().ToString("o")!=parentBirth || parent.HasExited) Environment.Exit(125);
                }}catch {Environment.Exit(125);}
                if(timer.Elapsed.TotalSeconds>45) Environment.Exit(124);
                System.Threading.Thread.Sleep(100);
            }
        });watchdog.IsBackground=true;watchdog.Start();
        try {return ObserveCore(loaderPath,transportPath,endpointPath,outputDirectory);}
        finally {System.Threading.Volatile.Write(ref done,true);watchdog.Join(1000);}
    }
    static object ObserveCore(string loaderPath, string transportPath, string endpointPath, string outputDirectory) {'''
s=one(s,anchor,wrapper)
s=one(s,'if (!amdInfo || !advertised.Contains("VK_KHR_shader_fma"))','if (!amdInfo || !advertised.Contains("VK_KHR_shader_fma") || !advertised.Contains("VK_KHR_shader_clock"))')
s=one(s,'            IntPtr fma = Allocate(32, owned); Marshal.WriteInt32(fma, 1000579000);','            IntPtr shaderClock = Allocate(24, owned); Marshal.WriteInt32(shaderClock,1000181000);\n            IntPtr fma = Allocate(32, owned); Marshal.WriteInt32(fma, 1000579000); Marshal.WriteIntPtr(fma,8,shaderClock);')
s=one(s,'            Marshal.WriteInt32(fma, 16, 0);','            uint subgroupClockSupported=U32(shaderClock,16),deviceClockSupported=U32(shaderClock,20);\n            if(subgroupClockSupported>1 || deviceClockSupported!=1 || U32(shaderClock,0)!=1000181000 || Marshal.ReadIntPtr(shaderClock,8)!=IntPtr.Zero) throw new Exception("Actual Device-clock feature unavailable or malformed");\n            Marshal.WriteInt32(shaderClock,16,0); Marshal.WriteInt32(shaderClock,20,1);\n            Marshal.WriteInt32(fma, 16, 0);')
s=one(s,'new string[] { "VK_KHR_shader_fma", "VK_AMD_shader_info" }','new string[] { "VK_KHR_shader_fma", "VK_AMD_shader_info", "VK_KHR_shader_clock" }')
s=one(s,'IntPtr enabledExtensions = Allocate(16, owned);','IntPtr enabledExtensions = Allocate(24, owned);')
s=one(s,'extensionCount = 2, extensions = enabledExtensions','extensionCount = 3, extensions = enabledExtensions')
s=one(s,'var setInfo = new SetLayoutInfo { sType = 32, count = 2','var setInfo = new SetLayoutInfo { sType = 32, count = 3')
s=one(s,'new Binding { binding = 1, type = 7, count = 1, stages = 0x20 } }, owned)','new Binding { binding = 1, type = 7, count = 1, stages = 0x20 }, new Binding { binding = 2, type = 7, count = 1, stages = 0x20 } }, owned)')
s=one(s,'b52e6373f6b30aac6c7fe877f6e07f3c47c7ac9eafa1671840eacd19960d7cb1','0eff0f81bf94582c5401e544b28bca9056ee77beb0264d4ece541aa8443dd076')
s=one(s,'785d94c53839693af8a0c752b8d7f0f1bd44789df17f3ab4eef79d3684faa559','1c4687976bb4efb00d2c33c55a9690d0128591ab22188903894e47217959576d')
s=one(s,'for (int i = 1; i < paths.Length; ++i)','for (int i = 0; i < paths.Length; ++i)')
s=one(s,'"pipelines_created", 1, "dispatches", 0','"logical_devices_created",1,"pipelines_created", 2, "dispatches", 0,"descriptor_bindings",new uint[]{0,1,2},"shader_device_clock_supported",deviceClockSupported,"shader_subgroup_clock_supported",subgroupClockSupported,"shader_device_clock_enabled",true,"shader_subgroup_clock_enabled",false')
s=one(s,'Endpoint-only unchanged production module/features/layout/flags with AMD tooling extension added. Static generated-code observations; no dispatch, timing comparison, numerical or full qualification verdict.','Exact unadopted native-default phase-observer modules; Device-clock feature and descriptorbinding2 explicitly enabled with AMD tooling extension. Empty pipeline cache/flags0 ordinary compilation; driver coldness unverified. Static generated-code query only; no dispatch, clock execution/rate, timing benefit, numerical/scientific/release verdict.')
save('native_clock_placement.cs',s)
child=(W/'prior-clock-child.ps1').read_text();child=one(child,'native_clock_features.cs','native_clock_placement.cs')
child=one(child," $result=[SiriusNativeClockFeatures]::Query('C:\\Windows\\System32\\vulkan-1.dll',[int]$owner.owner_pid,[string]$owner.owner_start_time_utc)"," $payloads=(New-Item -ItemType Directory -Path (Join-Path $Stage 'query-payloads')).FullName\n $result=[NativeClockPlacement]::Observe('C:\\Windows\\System32\\vulkan-1.dll',(Join-Path $Stage 'retained_transport.spv'),(Join-Path $Stage 'retained_endpoint.spv'),$payloads,[int]$owner.owner_pid,[string]$owner.owner_start_time_utc)")
child=child.replace('clock-features.json','clock-placement.json').replace('Native clock feature query completed','Native clock placement query completed');save('placement-child.ps1',child)
owner=(W/'prior-clock-owner.ps1').read_text();owner=owner.replace('clock-child.ps1','placement-child.ps1').replace('clock-features.json','clock-placement.json').replace('[int]$Seconds=45','[int]$Seconds=50').replace('6be3132b7023f9d6497ab8f67dc779cf8eb4a80c',head)
owner=one(owner,'$result.logical_devices_created -ne 0','$result.logical_devices_created -ne 1 -or $result.pipelines_created -ne 2 -or !$result.shader_device_clock_enabled')
owner=owner.replace('Read-only exact native shader-clock feature query; no logical device, pipeline, dispatch or scientific/release qualification. Fixed 45s/512MiB diagnostic guard.','Ordinary generated-code query for two exact diagnostic modules with Device-clock enabled and binding2; no dispatch/scientific/release qualification.50s native polling/sample512MiB childRSS guard;120s bridge plus separate cleanup bounds.')
save('placement-owner.ps1',owner)
save('emergency_native_cleanup.ps1',(W/'prior-emergency_native_cleanup.ps1').read_text())
bridge=(W/'prior-run_clock_query.py').read_text();bridge=bridge.replace("action='clock-query'","action='placement-query'").replace("WORK/'clock-query-bridge'","WORK/'placement-query-bridge'").replace('clock-owner.ps1','placement-owner.ps1').replace('clock-child.ps1','placement-child.ps1').replace('native_clock_features.cs','native_clock_placement.cs').replace('clock-features.json','clock-placement.json').replace("WORK/'prepare_clock_owner.py',WORK/'prepare_clock_bridge.py',WORK/'owner-preparation.json'","WORK/'prepare_placement.py',WORK/'preparation.json'").replace('native-shader-clock-features','native-shader-clock-placement')
bridge=bridge.replace("'-Seconds','45'","'-Seconds','50'").replace("owner['outer_guard_seconds']==45","owner['outer_guard_seconds']==50")
bridge=one(bridge,"assert result['logical_devices_created']==result['pipelines_created']==result['dispatches']==0","assert result['logical_devices_created']==1 and result['pipelines_created']==2 and result['dispatches']==0\n  assert result['shader_device_clock_supported']==1 and result['shader_device_clock_enabled'] and not result['shader_subgroup_clock_enabled'] and result['descriptor_bindings']==[0,1,2]\n  assert len(result['observations'])==2\n  for item in result['observations']:\n   stem='retained_'+item['stage'].lower();expected=seal(WORK/(stem+'.spv'))\n   assert item['module_bytes']==expected['bytes'] and item['module_sha256']==expected['sha256'] and item['pipeline_flags']==item['stage_flags']==item['pipeline_creation_result']==0\n   for q in item['queries']:\n    if q.get('captured'):\n     assert seal(WORK/'query-payloads'/q['file'])=={'bytes':q['returned_bytes'],'sha256':q['sha256']}\n    elif q.get('available'):assert q['requested_bytes']>8388608 and not q.get('captured',False)")
bridge=one(bridge,"assert result['extension_count']==len(result['extensions'])","assert result['enabled_extensions']==['VK_KHR_shader_fma','VK_AMD_shader_info','VK_KHR_shader_clock']")
bridge=bridge.replace("record['clock_features']=result['clock_features']","record['clock_placement']={'device_clock_enabled':True,'pipeline_count':2,'dispatches':0}").replace("'clock_features')","'clock_placement')")
save('run_placement_query.py',bridge);ast.parse(bridge)
inputs={str(p.relative_to(R)):seal(p) for p in [Path(__file__),W/'prior-native_endpoint_query.cs',W/'prior-clock-owner.ps1',W/'prior-clock-child.ps1',W/'prior-run_clock_query.py']}
save('preparation.json',json.dumps({'input_seals':inputs,'outputs':{n:seal(W/n) for n in ('native_clock_placement.cs','placement-owner.ps1','placement-child.ps1','run_placement_query.py','emergency_native_cleanup.ps1')},'scope':'Source preparation only. No helper/compiler/native invocation.'},indent=2)+'\n')
print(json.dumps({'outputs':json.loads((W/'preparation.json').read_text())['outputs']}))
