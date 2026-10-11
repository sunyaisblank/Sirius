using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Threading;

// Read-only Win64 feature query. No logical device or submitted GPU work.
public static class SiriusNativeClockFeatures {
    const string LoaderHash = "45ee2efc5c6d3986da8181775610f308122c752e55e703340d9733bbbc4b9f1a";
    const string IcdHash = "5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650";
    const string IcdPath = @"C:\WINDOWS\System32\DriverStore\FileRepository\amdvlk.inf_amd64_914ba89eaaafdf60\amdvlk64.dll";
    [StructLayout(LayoutKind.Sequential)] struct ApplicationInfo {
        public uint sType; public IntPtr pNext, applicationName; public uint applicationVersion;
        public IntPtr engineName; public uint engineVersion, apiVersion;
    }
    [StructLayout(LayoutKind.Sequential)] struct InstanceInfo {
        public uint sType; public IntPtr pNext; public uint flags; public IntPtr application;
        public uint layerCount; public IntPtr layers; public uint extensionCount; public IntPtr extensions;
    }
    [StructLayout(LayoutKind.Sequential)] struct ClockFeatures {
        public uint sType; public IntPtr pNext; public uint subgroupClock, deviceClock;
    }
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate IntPtr ProcAddress(IntPtr handle, [MarshalAs(UnmanagedType.LPStr)] string name);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreateInstance(ref InstanceInfo info, IntPtr allocator, out IntPtr instance);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void DestroyInstance(IntPtr instance, IntPtr allocator);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int EnumerateDevices(IntPtr instance, ref uint count, IntPtr devices);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void GetValues(IntPtr physical, IntPtr data);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int Extensions(IntPtr physical, IntPtr layer, ref uint count, IntPtr properties);
    [DllImport("kernel32.dll", CharSet=CharSet.Unicode, SetLastError=true)] static extern IntPtr LoadLibraryEx(string path, IntPtr file, uint flags);
    [DllImport("kernel32.dll", CharSet=CharSet.Ansi)] static extern IntPtr GetProcAddress(IntPtr module, string name);
    [DllImport("kernel32.dll")] static extern bool FreeLibrary(IntPtr module);
    static T Bind<T>(IntPtr address, string name) where T:class {
        if (address==IntPtr.Zero) throw new Exception("Missing Vulkan entry point: "+name);
        return (T)(object)Marshal.GetDelegateForFunctionPointer(address, typeof(T));
    }
    static Dictionary<string,object> Record(params object[] pairs) {
        var value=new Dictionary<string,object>();
        for(int i=0;i<pairs.Length;i+=2)value.Add((string)pairs[i],pairs[i+1]);return value;
    }
    static string Hash(byte[] bytes) {
        using(var sha=SHA256.Create())return BitConverter.ToString(sha.ComputeHash(bytes)).Replace("-","").ToLowerInvariant();
    }
    static uint U32(IntPtr p,int offset){return unchecked((uint)Marshal.ReadInt32(p,offset));}
    static IntPtr Allocate(int bytes,List<IntPtr> owned){
        if(bytes<=0 || bytes>1024*1024)throw new Exception("Allocation bound exceeded");
        var p=Marshal.AllocHGlobal(bytes);owned.Add(p);Marshal.Copy(new byte[bytes],0,p,bytes);return p;
    }
    static IntPtr Structure<T>(T value,List<IntPtr> owned) where T:struct {
        var p=Allocate(Marshal.SizeOf(typeof(T)),owned);Marshal.StructureToPtr(value,p,false);return p;
    }
    static void Check(int result,string operation){if(result!=0)throw new Exception(operation+" VkResult "+result);}
    static byte[] Bytes(IntPtr p,int count){var b=new byte[count];Marshal.Copy(p,b,0,count);return b;}
    public static object Query(string loaderPath,int parentPid,string parentBirth) {
        if(IntPtr.Size!=8 || !BitConverter.IsLittleEndian || Marshal.SizeOf(typeof(ApplicationInfo))!=48 || Marshal.SizeOf(typeof(InstanceInfo))!=64 || Marshal.SizeOf(typeof(ClockFeatures))!=24 || Marshal.OffsetOf(typeof(ClockFeatures),"subgroupClock").ToInt32()!=16 || Marshal.OffsetOf(typeof(ClockFeatures),"deviceClock").ToInt32()!=20)throw new Exception("Win64 ABI mismatch");
        if(Hash(File.ReadAllBytes(loaderPath))!=LoaderHash)throw new Exception("Pinned loader differs");
        bool done=false;var timer=Stopwatch.StartNew();
        var watchdog=new Thread(delegate(){
            while(!Volatile.Read(ref done)){
                try {using(var parent=Process.GetProcessById(parentPid)){
                    if(parent.StartTime.ToUniversalTime().ToString("o")!=parentBirth || parent.HasExited)Environment.Exit(125);
                }}catch {Environment.Exit(125);}
                if(timer.Elapsed.TotalSeconds>40)Environment.Exit(124);Thread.Sleep(100);
            }
        });watchdog.IsBackground=true;watchdog.Start();
        IntPtr loader=IntPtr.Zero,instance=IntPtr.Zero;DestroyInstance destroy=null;var owned=new List<IntPtr>();
        try {
            loader=LoadLibraryEx(loaderPath,IntPtr.Zero,8);if(loader==IntPtr.Zero)throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error());
            var proc=Bind<ProcAddress>(GetProcAddress(loader,"vkGetInstanceProcAddr"),"vkGetInstanceProcAddr");
            var create=Bind<CreateInstance>(proc(IntPtr.Zero,"vkCreateInstance"),"vkCreateInstance");
            var app=new ApplicationInfo{sType=0,apiVersion=(1u<<22)|(3u<<12)};
            var info=new InstanceInfo{sType=1,application=Structure(app,owned)};
            Check(create(ref info,IntPtr.Zero,out instance),"vkCreateInstance");
            destroy=Bind<DestroyInstance>(proc(instance,"vkDestroyInstance"),"vkDestroyInstance");
            var enumerate=Bind<EnumerateDevices>(proc(instance,"vkEnumeratePhysicalDevices"),"vkEnumeratePhysicalDevices");
            uint count=0;Check(enumerate(instance,ref count,IntPtr.Zero),"physical device count");if(count!=1)throw new Exception("Expected one native physical device");
            var devices=Allocate(8,owned);Check(enumerate(instance,ref count,devices),"complete physical device list");if(count!=1)throw new Exception("Physical device list changed");
            IntPtr physical=Marshal.ReadIntPtr(devices);
            var driver=Allocate(536,owned);Marshal.WriteInt32(driver,1000196000);
            var properties=Allocate(2048,owned);Marshal.WriteInt32(properties,1000059001);Marshal.WriteIntPtr(properties,8,driver);
            Bind<GetValues>(proc(instance,"vkGetPhysicalDeviceProperties2"),"vkGetPhysicalDeviceProperties2")(physical,properties);
            var basic=IntPtr.Add(properties,16);string deviceName=Marshal.PtrToStringAnsi(IntPtr.Add(basic,20));string driverInfo=Marshal.PtrToStringAnsi(IntPtr.Add(driver,276));
            if(deviceName!="AMD Radeon 780M Graphics" || driverInfo!="26.8.1 (LLPC)" || U32(basic,8)!=4098 || U32(basic,12)!=6400)throw new Exception("Native provider identity differs");
            var extensions=Bind<Extensions>(proc(instance,"vkEnumerateDeviceExtensionProperties"),"vkEnumerateDeviceExtensionProperties");
            uint extensionCount=0;Check(extensions(physical,IntPtr.Zero,ref extensionCount,IntPtr.Zero),"extension count");if(extensionCount==0 || extensionCount>1024)throw new Exception("Extension count bound");
            uint capacity=extensionCount;var extensionData=Allocate(checked((int)capacity*260),owned);Check(extensions(physical,IntPtr.Zero,ref extensionCount,extensionData),"complete extension list");if(capacity!=extensionCount)throw new Exception("Extension list changed");
            var names=new List<object>();bool advertised=false;uint version=0;
            for(uint i=0;i<extensionCount;i++){
                var item=IntPtr.Add(extensionData,checked((int)i*260));int length=0;while(length<256 && Marshal.ReadByte(item,length)!=0)length++;
                if(length==0 || length==256)throw new Exception("Malformed extension name");string name=Marshal.PtrToStringAnsi(item,length);uint v=U32(item,256);names.Add(Record("name",name,"spec_version",v));
                if(name=="VK_KHR_shader_clock"){if(advertised)throw new Exception("Duplicate clock extension");advertised=true;version=v;}
            }
            object features=null;
            if(advertised){
                // VkPhysicalDeviceFeatures2: header16 +55 VkBool32 +4 tail padding.
                var clock=Allocate(24,owned);Marshal.WriteInt32(clock,1000181000);
                var data=Allocate(240,owned);Marshal.WriteInt32(data,1000059000);Marshal.WriteIntPtr(data,8,clock);
                Bind<GetValues>(proc(instance,"vkGetPhysicalDeviceFeatures2"),"vkGetPhysicalDeviceFeatures2")(physical,data);
                uint subgroup=U32(clock,16),device=U32(clock,20);if(subgroup>1 || device>1)throw new Exception("Invalid VkBool32 feature");
                if(U32(clock,0)!=1000181000 || Marshal.ReadIntPtr(clock,8)!=IntPtr.Zero || U32(data,0)!=1000059000 || Marshal.ReadIntPtr(data,8)!=clock)throw new Exception("Feature chain headers changed");
                features=Record("shaderSubgroupClock",subgroup==1,"shaderDeviceClock",device==1,"shaderSubgroupClock_raw",subgroup,"shaderDeviceClock_raw",device,"clock_struct_bytes",Bytes(clock,24),"features2_struct_bytes",Bytes(data,240),"shaderInt64_supported",U32(data,16+160)!=0,"shaderFloat64_supported",U32(data,16+156)!=0,"enabled_features",0);
            }
            var modules=new List<object>();bool loaderLoaded=false,icdLoaded=false;
            foreach(ProcessModule module in Process.GetCurrentProcess().Modules){
                string name=Path.GetFileName(module.FileName);
                if(!name.Equals("vulkan-1.dll",StringComparison.OrdinalIgnoreCase) && !name.Equals("amdvlk64.dll",StringComparison.OrdinalIgnoreCase))continue;
                string hash=Hash(File.ReadAllBytes(module.FileName));
                if(name.Equals("vulkan-1.dll",StringComparison.OrdinalIgnoreCase)){
                    if(!Path.GetFullPath(module.FileName).Equals(Path.GetFullPath(loaderPath),StringComparison.OrdinalIgnoreCase) || hash!=LoaderHash)throw new Exception("Actual Vulkan loader differs");loaderLoaded=true;
                }else{
                    if(!module.FileName.Equals(IcdPath,StringComparison.OrdinalIgnoreCase) || hash!=IcdHash)throw new Exception("Actual AMD ICD differs");icdLoaded=true;
                }
                modules.Add(Record("name",name,"path",module.FileName,"bytes",new FileInfo(module.FileName).Length,"sha256",hash));
            }
            if(!loaderLoaded || !icdLoaded || modules.Count!=2)throw new Exception("Actual provider module census incomplete");
            return Record("completed",true,"instance_api","1.3","device",deviceName,"driver",driverInfo,"extensions",names,"extension_count",extensionCount,"clock_extension_advertised",advertised,"clock_extension_revision",version,"clock_features",features,"loaded_modules",modules,"logical_devices_created",0,"pipelines_created",0,"dispatches",0,"scope","Actual read-only feature admission only. No enabled logical device, compiled-placement, clock execution, timing, scientific, performance or release verdict.");
        }finally{
            try {if(instance!=IntPtr.Zero && destroy!=null)destroy(instance,IntPtr.Zero);}
            finally {foreach(var p in owned)Marshal.FreeHGlobal(p);if(loader!=IntPtr.Zero)FreeLibrary(loader);Volatile.Write(ref done,true);watchdog.Join(1000);}
        }
    }
}
