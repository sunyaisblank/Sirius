using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Security.Cryptography;

public static class SiriusShaderWave32Query {
  [DllImport("kernel32.dll", CharSet=CharSet.Unicode, SetLastError=true)] static extern IntPtr LoadLibraryEx(string path, IntPtr file, uint flags);
  [DllImport("kernel32.dll", CharSet=CharSet.Ansi)] static extern IntPtr GetProcAddress(IntPtr library, string name);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreateInstance(IntPtr info, IntPtr allocator, out IntPtr instance);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int EnumeratePhysical(IntPtr instance, ref uint count, IntPtr devices);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void PhysicalProperties(IntPtr physical, IntPtr properties);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void QueueProperties(IntPtr physical, ref uint count, IntPtr properties);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int EnumerateExtensions(IntPtr physical, IntPtr layer, ref uint count, IntPtr extensions);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreateDevice(IntPtr physical, IntPtr info, IntPtr allocator, out IntPtr device);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreateObject(IntPtr device, IntPtr info, IntPtr allocator, out ulong handle);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreatePipeline(IntPtr device, ulong cache, uint count, IntPtr info, IntPtr allocator, out ulong pipeline);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void DestroyObject(IntPtr device, ulong handle, IntPtr allocator);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void DestroyDispatchable(IntPtr handle, IntPtr allocator);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate IntPtr DeviceProc(IntPtr device, [MarshalAs(UnmanagedType.LPStr)] string name);
  [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int ShaderInfo(IntPtr device, ulong pipeline, uint stage, uint kind, ref UIntPtr bytes, IntPtr data);

  [StructLayout(LayoutKind.Sequential)] struct ApplicationInfo { public uint type; public IntPtr next, name; public uint version; public IntPtr engine; public uint engineVersion, api; }
  [StructLayout(LayoutKind.Sequential)] struct InstanceInfo { public uint type; public IntPtr next; public uint flags; public IntPtr application; public uint layerCount; public IntPtr layers; public uint extensionCount; public IntPtr extensions; }
  [StructLayout(LayoutKind.Sequential)] struct QueueInfo { public uint type; public IntPtr next; public uint flags, family, count; public IntPtr priorities; }
  [StructLayout(LayoutKind.Sequential)] struct DeviceInfo { public uint type; public IntPtr next; public uint flags, queueCount; public IntPtr queues; public uint layerCount; public IntPtr layers; public uint extensionCount; public IntPtr extensions, features; }
  [StructLayout(LayoutKind.Sequential)] struct ModuleInfo { public uint type; public IntPtr next; public uint flags; public UIntPtr bytes; public IntPtr code; }
  [StructLayout(LayoutKind.Sequential)] struct Binding { public uint binding, descriptorType, count, stages; public IntPtr samplers; }
  [StructLayout(LayoutKind.Sequential)] struct SetInfo { public uint type; public IntPtr next; public uint flags, bindingCount; public IntPtr bindings; }
  [StructLayout(LayoutKind.Sequential)] struct LayoutInfo { public uint type; public IntPtr next; public uint flags, setCount; public IntPtr sets; public uint pushCount; public IntPtr pushes; }
  [StructLayout(LayoutKind.Sequential)] struct StageInfo { public uint type; public IntPtr next; public uint flags, stage; public ulong module; public IntPtr name, specialization; }
  [StructLayout(LayoutKind.Sequential)] struct PipelineInfo { public uint type; public IntPtr next; public uint flags; public StageInfo stage; public ulong layout, basePipeline; public int baseIndex; }
  [StructLayout(LayoutKind.Sequential)] struct ResourceUsage { public uint vgprs, sgprs, ldsPerGroup; public UIntPtr ldsBytes, scratchBytes; }
  [StructLayout(LayoutKind.Sequential)] struct Statistics { public uint stages; public ResourceUsage resource; public uint physicalVgprs, physicalSgprs, availableVgprs, availableSgprs, groupX, groupY, groupZ; }

  static void Need(bool condition, string message) { if (!condition) throw new InvalidOperationException(message); }
  static void Ok(int result, string operation) { Need(result == 0, operation + " VkResult=" + result); }
  static T Function<T>(IntPtr library, string name) where T:class { IntPtr p=GetProcAddress(library,name); Need(p!=IntPtr.Zero,"Missing export "+name); return (T)(object)Marshal.GetDelegateForFunctionPointer(p,typeof(T)); }
  static string Hash(byte[] bytes) { using (SHA256 h=SHA256.Create()) return BitConverter.ToString(h.ComputeHash(bytes)).Replace("-","").ToLowerInvariant(); }
  static IntPtr Allocate(List<IntPtr> memory, int bytes) { IntPtr p=Marshal.AllocHGlobal(bytes); memory.Add(p); Marshal.Copy(new byte[bytes],0,p,bytes); return p; }
  static IntPtr Structure<T>(List<IntPtr> memory, T value) where T:struct { IntPtr p=Allocate(memory,Marshal.SizeOf(typeof(T))); Marshal.StructureToPtr(value,p,false); return p; }
  static IntPtr Text(List<IntPtr> memory, string value) { IntPtr p=Marshal.StringToHGlobalAnsi(value); memory.Add(p); return p; }
  static Dictionary<string,object> Fields(params object[] pairs) { var d=new Dictionary<string,object>(); for(int i=0;i<pairs.Length;i+=2)d.Add((string)pairs[i],pairs[i+1]); return d; }

  public static Dictionary<string,object> Run(string loader, string shader, string output) {
    Need(IntPtr.Size==8,"x64 required");
    Need(Marshal.SizeOf(typeof(ApplicationInfo))==48 && Marshal.SizeOf(typeof(InstanceInfo))==64 && Marshal.SizeOf(typeof(QueueInfo))==40 && Marshal.SizeOf(typeof(DeviceInfo))==72,"Instance/device ABI mismatch");
    Need(Marshal.SizeOf(typeof(ModuleInfo))==40 && Marshal.SizeOf(typeof(Binding))==24 && Marshal.SizeOf(typeof(SetInfo))==32 && Marshal.SizeOf(typeof(LayoutInfo))==48,"Shader/layout ABI mismatch");
    Need(Marshal.SizeOf(typeof(StageInfo))==48 && Marshal.SizeOf(typeof(PipelineInfo))==96 && Marshal.SizeOf(typeof(Statistics))==72,"Pipeline/statistics ABI mismatch");
    var memory=new List<IntPtr>(); IntPtr instance=IntPtr.Zero, device=IntPtr.Zero;
    ulong module=0, set=0, layout=0, pipeline=0;
    IntPtr library=LoadLibraryEx(loader,IntPtr.Zero,0x1100); Need(library!=IntPtr.Zero,"SDK loader load failed "+Marshal.GetLastWin32Error());
    try {
      var app=new ApplicationInfo { type=0, name=Text(memory,"Sirius shader introspection"), engine=Text(memory,"Sirius"), api=(1u<<22)|(2u<<12) };
      var info=new InstanceInfo { type=1, application=Structure(memory,app) };
      Ok(Function<CreateInstance>(library,"vkCreateInstance")(Structure(memory,info),IntPtr.Zero,out instance),"CreateInstance");
      uint count=0; var enumerate=Function<EnumeratePhysical>(library,"vkEnumeratePhysicalDevices");
      Ok(enumerate(instance,ref count,IntPtr.Zero),"Physical count"); Need(count>0 && count<=32,"Physical count bound");
      IntPtr physicals=Allocate(memory,checked((int)count*8)); Ok(enumerate(instance,ref count,physicals),"Physical devices");
      IntPtr physical=IntPtr.Zero; uint driver=0; var names=new List<string>();
      for(int i=0;i<count;i++) { IntPtr candidate=Marshal.ReadIntPtr(physicals,i*8), properties=Allocate(memory,4096); Function<PhysicalProperties>(library,"vkGetPhysicalDeviceProperties")(candidate,properties); string name=Marshal.PtrToStringAnsi(IntPtr.Add(properties,20)); names.Add(name); if(name=="AMD Radeon 780M Graphics") { Need(physical==IntPtr.Zero,"Ambiguous AMD device"); physical=candidate; driver=(uint)Marshal.ReadInt32(properties,4); } }
      Need(physical!=IntPtr.Zero,"Required native device absent");
      uint extensions=0; var enumerateExtensions=Function<EnumerateExtensions>(library,"vkEnumerateDeviceExtensionProperties");
      Ok(enumerateExtensions(physical,IntPtr.Zero,ref extensions,IntPtr.Zero),"Extension count"); Need(extensions<=4096,"Extension count bound");
      IntPtr extensionRows=Allocate(memory,checked((int)extensions*260)); Ok(enumerateExtensions(physical,IntPtr.Zero,ref extensions,extensionRows),"Extensions"); bool shaderInfo=false, subgroupControl=false;
      for(int i=0;i<extensions;i++) {
        string extension=Marshal.PtrToStringAnsi(IntPtr.Add(extensionRows,i*260));
        if(extension=="VK_AMD_shader_info")shaderInfo=true;
        if(extension=="VK_EXT_subgroup_size_control")subgroupControl=true;
      }
      Need(subgroupControl,"VK_EXT_subgroup_size_control unavailable");
      IntPtr subgroupFeatures=Allocate(memory,24); Marshal.WriteInt32(subgroupFeatures,0,1000225002);
      IntPtr features2=Allocate(memory,512); Marshal.WriteInt32(features2,0,1000059000); Marshal.WriteIntPtr(features2,8,subgroupFeatures);
      Function<PhysicalProperties>(library,"vkGetPhysicalDeviceFeatures2")(physical,features2);
      Need(Marshal.ReadInt32(subgroupFeatures,16)!=0,"subgroupSizeControl feature absent");
      // Enable only subgroupSizeControl, without computeFullSubgroups or other optional features.
      Marshal.WriteInt32(subgroupFeatures,16,1); Marshal.WriteInt32(subgroupFeatures,20,0);
      IntPtr subgroupProperties=Allocate(memory,32); Marshal.WriteInt32(subgroupProperties,0,1000225000);
      IntPtr properties2=Allocate(memory,4096); Marshal.WriteInt32(properties2,0,1000059001); Marshal.WriteIntPtr(properties2,8,subgroupProperties);
      Function<PhysicalProperties>(library,"vkGetPhysicalDeviceProperties2")(physical,properties2);
      uint minimum=(uint)Marshal.ReadInt32(subgroupProperties,16), maximum=(uint)Marshal.ReadInt32(subgroupProperties,20);
      uint maximumSubgroups=(uint)Marshal.ReadInt32(subgroupProperties,24), requiredStages=(uint)Marshal.ReadInt32(subgroupProperties,28);
      Need(minimum<=32 && maximum>=32 && maximumSubgroups>=2 && (requiredStages&32)!=0,"Required32 subgroup compute support absent");
      Need(shaderInfo,"VK_AMD_shader_info unavailable");
      uint families=0; var queues=Function<QueueProperties>(library,"vkGetPhysicalDeviceQueueFamilyProperties"); queues(physical,ref families,IntPtr.Zero); Need(families>0 && families<=128,"Queue count bound");
      IntPtr familyRows=Allocate(memory,checked((int)families*24)); queues(physical,ref families,familyRows); int family=-1;
      for(int i=0;i<families;i++)if((Marshal.ReadInt32(familyRows,i*24)&2)!=0 && Marshal.ReadInt32(familyRows,i*24+4)>0) { family=i;break; }
      Need(family>=0,"Compute queue absent");
      IntPtr priority=Allocate(memory,4); Marshal.Copy(new float[]{1.0f},0,priority,1);
      var queue=new QueueInfo { type=2, family=(uint)family, count=1, priorities=priority };
      IntPtr extensionPointer=Allocate(memory,16); Marshal.WriteIntPtr(extensionPointer,0,Text(memory,"VK_AMD_shader_info")); Marshal.WriteIntPtr(extensionPointer,8,Text(memory,"VK_EXT_subgroup_size_control"));
      var deviceInfo=new DeviceInfo { type=3, next=subgroupFeatures, queueCount=1, queues=Structure(memory,queue), extensionCount=2, extensions=extensionPointer };
      Ok(Function<CreateDevice>(library,"vkCreateDevice")(physical,Structure(memory,deviceInfo),IntPtr.Zero,out device),"CreateDevice");
      IntPtr queryPointer=Function<DeviceProc>(library,"vkGetDeviceProcAddr")(device,"vkGetShaderInfoAMD"); Need(queryPointer!=IntPtr.Zero,"Shader query pointer absent");
      var query=(ShaderInfo)Marshal.GetDelegateForFunctionPointer(queryPointer,typeof(ShaderInfo));
      IntPtr bindings=Allocate(memory,48); for(int i=0;i<2;i++)Marshal.StructureToPtr(new Binding { binding=(uint)i, descriptorType=7, count=1, stages=32 },IntPtr.Add(bindings,i*24),false);
      var setInfo=new SetInfo { type=32, bindingCount=2, bindings=bindings };
      Ok(Function<CreateObject>(library,"vkCreateDescriptorSetLayout")(device,Structure(memory,setInfo),IntPtr.Zero,out set),"Descriptor layout");
      IntPtr sets=Allocate(memory,8); Marshal.WriteInt64(sets,unchecked((long)set)); var layoutInfo=new LayoutInfo { type=30, setCount=1, sets=sets };
      Ok(Function<CreateObject>(library,"vkCreatePipelineLayout")(device,Structure(memory,layoutInfo),IntPtr.Zero,out layout),"Pipeline layout");
      byte[] code=File.ReadAllBytes(shader); Need(code.Length>20 && code.Length<=1048576 && code.Length%4==0 && BitConverter.ToUInt32(code,0)==0x07230203,"SPIR-V input bound/magic");
      IntPtr codePointer=Allocate(memory,code.Length); Marshal.Copy(code,0,codePointer,code.Length); var moduleInfo=new ModuleInfo { type=16, bytes=new UIntPtr((uint)code.Length), code=codePointer };
      Ok(Function<CreateObject>(library,"vkCreateShaderModule")(device,Structure(memory,moduleInfo),IntPtr.Zero,out module),"Shader module");
      IntPtr requiredSize=Allocate(memory,24); Marshal.WriteInt32(requiredSize,0,1000225001); Marshal.WriteInt32(requiredSize,16,32);
      var pipelineInfo=new PipelineInfo { type=29, stage=new StageInfo { type=18, next=requiredSize, stage=32, module=module, name=Text(memory,"main") }, layout=layout, baseIndex=0 };
      Ok(Function<CreatePipeline>(library,"vkCreateComputePipelines")(device,0,1,Structure(memory,pipelineInfo),IntPtr.Zero,out pipeline),"Compute pipeline");
      Directory.CreateDirectory(output); var outputs=new List<object>();
      for(uint kind=0;kind<3;kind++) { UIntPtr bytes=UIntPtr.Zero; int result=query(device,pipeline,32,kind,ref bytes,IntPtr.Zero); ulong size=bytes.ToUInt64(); var row=Fields("kind",kind,"size_result",result,"requested_bytes",size); outputs.Add(row); if(result!=0)continue; Need(size>0 && size<=16777216,"Shader query output bound"); IntPtr data=Allocate(memory,checked((int)size)); result=query(device,pipeline,32,kind,ref bytes,data); row.Add("data_result",result); row.Add("returned_bytes",bytes.ToUInt64()); if(result!=0)continue; Need(bytes.ToUInt64()<=size,"Shader query grew beyond allocation"); byte[] payload=new byte[checked((int)bytes.ToUInt64())]; Marshal.Copy(data,payload,0,payload.Length); string path=Path.Combine(output,"shader-info-"+kind+(kind==2?".isa.txt":".bin")); File.WriteAllBytes(path,payload); row.Add("path",path); row.Add("sha256",Hash(payload)); if(kind==0) { Need(payload.Length==72,"Unexpected statistics ABI"); var s=(Statistics)Marshal.PtrToStructure(data,typeof(Statistics)); row.Add("statistics",Fields("stage_mask",s.stages,"used_vgprs",s.resource.vgprs,"used_sgprs",s.resource.sgprs,"lds_per_group",s.resource.ldsPerGroup,"lds_bytes",s.resource.ldsBytes.ToUInt64(),"scratch_bytes",s.resource.scratchBytes.ToUInt64(),"physical_vgprs",s.physicalVgprs,"physical_sgprs",s.physicalSgprs,"available_vgprs",s.availableVgprs,"available_sgprs",s.availableSgprs,"workgroup",new uint[]{s.groupX,s.groupY,s.groupZ})); } }
      var modules=new List<object>(); foreach(ProcessModule loaded in Process.GetCurrentProcess().Modules)if(loaded.ModuleName.Equals("vulkan-1.dll",StringComparison.OrdinalIgnoreCase)||loaded.ModuleName.Equals("amdvlk64.dll",StringComparison.OrdinalIgnoreCase)) { byte[] b=File.ReadAllBytes(loaded.FileName); modules.Add(Fields("name",loaded.ModuleName,"path",loaded.FileName,"bytes",b.Length,"sha256",Hash(b),"version",loaded.FileVersionInfo.FileVersion)); }
      return Fields("device","AMD Radeon 780M Graphics","enumerated_devices",names,"driver_version_raw",driver,"queue_family",family,"loader_requested",loader,"loaded_modules",modules,"shader",shader,"shader_bytes",code.Length,"shader_sha256",Hash(code),"entry_point","main","storage_bindings",2,"pipeline_flags",0,"pipeline_cache","VK_NULL_HANDLE","enabled_device_extensions",new string[]{"VK_AMD_shader_info","VK_EXT_subgroup_size_control"},"optional_device_features_enabled",new string[]{"subgroupSizeControl"},"required_subgroup_size",32,"queried_subgroup_properties",Fields("minimum",minimum,"maximum",maximum,"maximum_subgroups",maximumSubgroups,"required_stages",requiredStages),"query_results",outputs,"dispatches",0,"numerical_or_performance_qualification",false);
    } finally {
      if(device!=IntPtr.Zero) {
        if(pipeline!=0)Function<DestroyObject>(library,"vkDestroyPipeline")(device,pipeline,IntPtr.Zero);
        if(module!=0)Function<DestroyObject>(library,"vkDestroyShaderModule")(device,module,IntPtr.Zero);
        if(layout!=0)Function<DestroyObject>(library,"vkDestroyPipelineLayout")(device,layout,IntPtr.Zero);
        if(set!=0)Function<DestroyObject>(library,"vkDestroyDescriptorSetLayout")(device,set,IntPtr.Zero);
        Function<DestroyDispatchable>(library,"vkDestroyDevice")(device,IntPtr.Zero);
      }
      if(instance!=IntPtr.Zero)Function<DestroyDispatchable>(library,"vkDestroyInstance")(instance,IntPtr.Zero);
      for(int i=memory.Count-1;i>=0;i--)Marshal.FreeHGlobal(memory[i]);
      // The process owns the SDK library reference until exit, after all Vulkan objects are destroyed.
    }
  }
}
