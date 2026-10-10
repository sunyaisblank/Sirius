// Finite stateful raw/reference gate for isolated native paired Transport.
// Reads exact validated modules; this is not complete product qualification.
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;

public static class NativePairedTransportGate {
    const string DllHash = "45ee2efc5c6d3986da8181775610f308122c752e55e703340d9733bbbc4b9f1a";
    const ulong AllocationLimit = 8 * 1024 * 1024;
    const ulong FenceTimeoutNanoseconds = 1000000000;

    [StructLayout(LayoutKind.Sequential)] struct ApplicationInfo {
        public uint sType; public IntPtr pNext, applicationName; public uint applicationVersion;
        public IntPtr engineName; public uint engineVersion, apiVersion;
    }
    [StructLayout(LayoutKind.Sequential)] struct InstanceInfo {
        public uint sType; public IntPtr pNext; public uint flags; public IntPtr application;
        public uint layerCount; public IntPtr layers; public uint extensionCount; public IntPtr extensions;
    }
    [StructLayout(LayoutKind.Sequential)] struct QueueInfo {
        public uint sType; public IntPtr pNext; public uint flags, family, count; public IntPtr priorities;
    }
    [StructLayout(LayoutKind.Sequential)] struct DeviceInfo {
        public uint sType; public IntPtr pNext; public uint flags, queueCount; public IntPtr queues;
        public uint layerCount; public IntPtr layers; public uint extensionCount;
        public IntPtr extensions, features;
    }
    [StructLayout(LayoutKind.Sequential)] struct BufferInfo {
        public uint sType; public IntPtr pNext; public uint flags; public ulong size;
        public uint usage, sharing, familyCount; public IntPtr families;
    }
    [StructLayout(LayoutKind.Sequential)] struct Requirements {
        public ulong size, alignment; public uint memoryTypeBits;
    }
    [StructLayout(LayoutKind.Sequential)] struct AllocationInfo {
        public uint sType; public IntPtr pNext; public ulong size; public uint memoryType;
    }
    [StructLayout(LayoutKind.Sequential)] struct ShaderInfo {
        public uint sType; public IntPtr pNext; public uint flags; public UIntPtr codeSize; public IntPtr code;
    }
    [StructLayout(LayoutKind.Sequential)] struct Binding {
        public uint binding, type, count, stages; public IntPtr immutableSamplers;
    }
    [StructLayout(LayoutKind.Sequential)] struct SetLayoutInfo {
        public uint sType; public IntPtr pNext; public uint flags, count; public IntPtr bindings;
    }
    [StructLayout(LayoutKind.Sequential)] struct PipelineLayoutInfo {
        public uint sType; public IntPtr pNext; public uint flags, setCount; public IntPtr sets;
        public uint pushCount; public IntPtr pushes;
    }
    [StructLayout(LayoutKind.Sequential)] struct ShaderStage {
        public uint sType; public IntPtr pNext; public uint flags, stage; public ulong module;
        public IntPtr name, specialization;
    }
    [StructLayout(LayoutKind.Sequential)] struct PipelineInfo {
        public uint sType; public IntPtr pNext; public uint flags; public ShaderStage stage;
        public ulong layout, basePipeline; public int baseIndex;
    }
    [StructLayout(LayoutKind.Sequential)] struct PoolSize { public uint type, count; }
    [StructLayout(LayoutKind.Sequential)] struct PoolInfo {
        public uint sType; public IntPtr pNext; public uint flags, maxSets, count; public IntPtr sizes;
    }
    [StructLayout(LayoutKind.Sequential)] struct SetAllocationInfo {
        public uint sType; public IntPtr pNext; public ulong pool; public uint count; public IntPtr layouts;
    }
    [StructLayout(LayoutKind.Sequential)] struct DescriptorBuffer { public ulong buffer, offset, range; }
    [StructLayout(LayoutKind.Sequential)] struct WriteSet {
        public uint sType; public IntPtr pNext; public ulong set; public uint binding, element, count, type;
        public IntPtr images, buffers, views;
    }
    [StructLayout(LayoutKind.Sequential)] struct CommandPoolInfo {
        public uint sType; public IntPtr pNext; public uint flags, family;
    }
    [StructLayout(LayoutKind.Sequential)] struct CommandAllocationInfo {
        public uint sType; public IntPtr pNext; public ulong pool; public uint level, count;
    }
    [StructLayout(LayoutKind.Sequential)] struct CommandBeginInfo {
        public uint sType; public IntPtr pNext; public uint flags; public IntPtr inheritance;
    }
    [StructLayout(LayoutKind.Sequential)] struct MemoryBarrier {
        public uint sType; public IntPtr pNext; public uint source, destination;
    }
    [StructLayout(LayoutKind.Sequential)] struct FenceInfo {
        public uint sType; public IntPtr pNext; public uint flags;
    }
    [StructLayout(LayoutKind.Sequential)] struct QueryPoolInfo {
        public uint sType; public IntPtr pNext; public uint flags, type, count, statistics;
    }
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void ResetQueries(IntPtr command, ulong pool, uint first, uint count);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void Timestamp(IntPtr command, uint stage, ulong pool, uint query);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int QueryResults(IntPtr device, ulong pool, uint first, uint count, UIntPtr bytes, IntPtr results, ulong stride, uint flags);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int ResetCommand(IntPtr command, uint flags);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int ResetFences(IntPtr device, uint count, IntPtr fences);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int QueueIdle(IntPtr queue);
    [StructLayout(LayoutKind.Sequential)] struct SubmitInfo {
        public uint sType; public IntPtr pNext; public uint waitCount; public IntPtr waits, waitStages;
        public uint commandCount; public IntPtr commands; public uint signalCount; public IntPtr signals;
    }

    // All dispatchable handles are IntPtr; Win64 non-dispatchable handles are ulong.
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate IntPtr ProcAddress(IntPtr handle, [MarshalAs(UnmanagedType.LPStr)] string name);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreateInstance(ref InstanceInfo info, IntPtr allocator, out IntPtr instance);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void DestroyInstance(IntPtr instance, IntPtr allocator);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int EnumerateDevices(IntPtr instance, ref uint count, IntPtr devices);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void Properties(IntPtr physical, IntPtr properties);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void Families(IntPtr physical, ref uint count, IntPtr properties);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreateDevice(IntPtr physical, ref DeviceInfo info, IntPtr allocator, out IntPtr device);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void DestroyDevice(IntPtr device, IntPtr allocator);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void GetQueue(IntPtr device, uint family, uint index, out IntPtr queue);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreateObject(IntPtr device, IntPtr info, IntPtr allocator, out ulong value);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void DestroyObject(IntPtr device, ulong value, IntPtr allocator);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void BufferRequirements(IntPtr device, ulong buffer, out Requirements requirements);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int BindMemory(IntPtr device, ulong buffer, ulong memory, ulong offset);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int MapMemory(IntPtr device, ulong memory, ulong offset, ulong size, uint flags, out IntPtr mapped);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void UnmapMemory(IntPtr device, ulong memory);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int CreatePipelines(IntPtr device, ulong cache, uint count, IntPtr infos, IntPtr allocator, IntPtr pipelines);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int AllocateSets(IntPtr device, IntPtr info, IntPtr sets);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void UpdateSets(IntPtr device, uint count, IntPtr writes, uint copyCount, IntPtr copies);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int AllocateCommands(IntPtr device, IntPtr info, out IntPtr command);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int BeginCommand(IntPtr command, ref CommandBeginInfo info);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int EndCommand(IntPtr command);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void Barrier(IntPtr command, uint source, uint destination, uint flags, uint memoryCount, ref MemoryBarrier memory, uint bufferCount, IntPtr buffers, uint imageCount, IntPtr images);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void BindPipeline(IntPtr command, uint point, ulong pipeline);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void BindSets(IntPtr command, uint point, ulong layout, uint first, uint count, IntPtr sets, uint offsetCount, IntPtr offsets);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate void Dispatch(IntPtr command, uint x, uint y, uint z);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int Submit(IntPtr queue, uint count, ref SubmitInfo info, ulong fence);
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int WaitFences(IntPtr device, uint count, IntPtr fences, uint waitAll, ulong timeout);
    [DllImport("kernel32.dll", CharSet=CharSet.Unicode, SetLastError=true)] static extern IntPtr LoadLibraryEx(string path, IntPtr file, uint flags);
    [DllImport("kernel32.dll", CharSet=CharSet.Ansi)] static extern IntPtr GetProcAddress(IntPtr module, string name);
    [DllImport("kernel32.dll")] static extern bool FreeLibrary(IntPtr module);
    [DllImport("kernel32.dll", CharSet=CharSet.Unicode)] static extern uint GetModuleFileName(IntPtr module, StringBuilder path, int count);

    static T Bind<T>(IntPtr address, string name) where T : class {
        if (address == IntPtr.Zero) throw new Exception("Missing entry point " + name);
        return (T)(object)Marshal.GetDelegateForFunctionPointer(address, typeof(T));
    }
    static Dictionary<string, object> Record(params object[] pairs) {
        var result = new Dictionary<string, object>();
        for (int i = 0; i < pairs.Length; i += 2) result.Add((string)pairs[i], pairs[i + 1]);
        return result;
    }
    static string Hash(byte[] bytes) {
        using (var sha = SHA256.Create()) return BitConverter.ToString(sha.ComputeHash(bytes)).Replace("-", "").ToLowerInvariant();
    }
    static byte[] ReadBounded(string path) {
        long length = new FileInfo(path).Length;
        if (length <= 0 || length > 8 * 1024 * 1024 || (length & 3) != 0) throw new Exception("Invalid bounded word file");
        byte[] bytes = File.ReadAllBytes(path);
        if (bytes.LongLength != length) throw new Exception("File changed during read");
        return bytes;
    }
    static uint Word(byte[] bytes, int offset) { return BitConverter.ToUInt32(bytes, offset); }
    static uint U32(IntPtr pointer, int offset) { return unchecked((uint)Marshal.ReadInt32(pointer, offset)); }
    static ulong U64(IntPtr pointer, int offset) { return unchecked((ulong)Marshal.ReadInt64(pointer, offset)); }
    static void Check(int result, string operation) {
        if (result != 0) throw new Exception(operation + " VkResult " + result);
    }
    static IntPtr Allocate(int size, List<IntPtr> owned) {
        IntPtr pointer = Marshal.AllocHGlobal(size);
        owned.Add(pointer); Marshal.Copy(new byte[size], 0, pointer, size); return pointer;
    }
    static IntPtr Structure<T>(T value, List<IntPtr> owned) where T : struct {
        IntPtr pointer = Allocate(Marshal.SizeOf(typeof(T)), owned);
        Marshal.StructureToPtr(value, pointer, false); return pointer;
    }
    static IntPtr Structures<T>(T[] values, List<IntPtr> owned) where T : struct {
        int stride = Marshal.SizeOf(typeof(T));
        IntPtr pointer = Allocate(checked(stride * values.Length), owned);
        for (int i = 0; i < values.Length; ++i) Marshal.StructureToPtr(values[i], IntPtr.Add(pointer, checked(stride * i)), false);
        return pointer;
    }
    static void Abi<T>(int expected) where T : struct {
        if (Marshal.SizeOf(typeof(T)) != expected) throw new Exception("Unexpected Win64 ABI: " + typeof(T).Name);
    }
    static void VerifyAbi() {
        if (IntPtr.Size != 8 || !BitConverter.IsLittleEndian) throw new Exception("Little-endian Win64 required");
        Abi<ApplicationInfo>(48); Abi<InstanceInfo>(64); Abi<QueueInfo>(40); Abi<DeviceInfo>(72);
        Abi<BufferInfo>(56); Abi<Requirements>(24); Abi<AllocationInfo>(32); Abi<ShaderInfo>(40);
        Abi<Binding>(24); Abi<SetLayoutInfo>(32); Abi<PipelineLayoutInfo>(48); Abi<ShaderStage>(48);
        Abi<PipelineInfo>(96); Abi<PoolSize>(8); Abi<PoolInfo>(40); Abi<SetAllocationInfo>(40);
        Abi<DescriptorBuffer>(24); Abi<WriteSet>(64); Abi<CommandPoolInfo>(24);
        Abi<CommandAllocationInfo>(32); Abi<CommandBeginInfo>(32); Abi<MemoryBarrier>(24);
        Abi<FenceInfo>(24); Abi<SubmitInfo>(72); Abi<QueryPoolInfo>(32);
    }

    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int DeviceExtensions(IntPtr physical,IntPtr layer,ref uint count,IntPtr properties);
    public static object Run(string loaderPath, string modulePath, string inputPath,
                             string expectedPath, string actualPath, uint Rows, int outputWords,
                             string moduleBinding, string inputBinding, bool wide, string alternativePath, string alternativeBinding,
                             string[] sequencePaths, string[] sequenceBindings, uint[] activeCounts,
                             int[] uploadBytes, string[] sequenceNames) {
        File.WriteAllText(actualPath + ".setup.txt", "");
        Action<string> checkpoint = delegate(string stage) {
            File.AppendAllText(actualPath + ".setup.txt", DateTime.UtcNow.Ticks + "|" +
                System.Diagnostics.Stopwatch.GetTimestamp() + "|" + stage + Environment.NewLine);
        };
        checkpoint("entered");
        VerifyAbi();
        if (!Path.IsPathRooted(loaderPath)) throw new Exception("Absolute pinned DLL path required");
        byte[] dll = ReadBounded(loaderPath), module = ReadBounded(modulePath);
        byte[] input = ReadBounded(inputPath);
        byte[] alternative = ReadBounded(alternativePath);
        if (Hash(alternative) != alternativeBinding || Word(alternative, 0) != 0x07230203)
            throw new Exception("Exact alternative Transport module binding mismatch");
        bool compare = !String.IsNullOrEmpty(expectedPath);
        byte[] expected = compare ? ReadBounded(expectedPath) : new byte[checked(outputWords * 4)];
        string moduleHash = Hash(module), inputHash = Hash(input);
        if (Hash(dll) != DllHash || inputHash != inputBinding || moduleHash != moduleBinding)
            throw new Exception("Exact loader/module/input binding mismatch");
        if (Word(module, 0) != 0x07230203 || Word(input, 0) != Rows ||
            Rows == 0 || Rows > 128 || expected.Length != checked(outputWords * 4))
            throw new Exception("Bounded Transport corpus shape mismatch");
        checkpoint("input-bindings-validated");

        IntPtr loader = LoadLibraryEx(loaderPath, IntPtr.Zero, 8);
        if (loader == IntPtr.Zero) throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error());
        var owned = new List<IntPtr>(); var cleanup = new List<Action>();
        bool submitted = false, completed = false;
        try {
            var loadedPath = new StringBuilder(4096);
            uint pathLength = GetModuleFileName(loader, loadedPath, loadedPath.Capacity);
            if (pathLength == 0 || pathLength >= loadedPath.Capacity ||
                !String.Equals(Path.GetFullPath(loadedPath.ToString()), Path.GetFullPath(loaderPath), StringComparison.OrdinalIgnoreCase) ||
                Hash(ReadBounded(loadedPath.ToString())) != DllHash)
                throw new Exception("Loaded DLL path/hash differs from requested pinned runtime");
            var instanceProc = Bind<ProcAddress>(GetProcAddress(loader, "vkGetInstanceProcAddr"), "vkGetInstanceProcAddr");
            var createInstance = Bind<CreateInstance>(instanceProc(IntPtr.Zero, "vkCreateInstance"), "vkCreateInstance");
            var app = new ApplicationInfo { sType = 0, apiVersion = (1u << 22) | (3u << 12) };
            var instanceInfo = new InstanceInfo { sType = 1, application = Structure(app, owned) };
            IntPtr instance; Check(createInstance(ref instanceInfo, IntPtr.Zero, out instance), "vkCreateInstance");
            var destroyInstance = Bind<DestroyInstance>(instanceProc(instance, "vkDestroyInstance"), "vkDestroyInstance");
            cleanup.Add(delegate { destroyInstance(instance, IntPtr.Zero); });
            var enumerate = Bind<EnumerateDevices>(instanceProc(instance, "vkEnumeratePhysicalDevices"), "vkEnumeratePhysicalDevices");
            var properties = Bind<Properties>(instanceProc(instance, "vkGetPhysicalDeviceProperties2"), "vkGetPhysicalDeviceProperties2");
            var memoryProperties = Bind<Properties>(instanceProc(instance, "vkGetPhysicalDeviceMemoryProperties"), "vkGetPhysicalDeviceMemoryProperties");
            var families = Bind<Families>(instanceProc(instance, "vkGetPhysicalDeviceQueueFamilyProperties"), "vkGetPhysicalDeviceQueueFamilyProperties");
            var createDevice = Bind<CreateDevice>(instanceProc(instance, "vkCreateDevice"), "vkCreateDevice");
            var deviceProc = Bind<ProcAddress>(instanceProc(instance, "vkGetDeviceProcAddr"), "vkGetDeviceProcAddr");
            uint physicalCount = 0; Check(enumerate(instance, ref physicalCount, IntPtr.Zero), "physical count");
            if (physicalCount != 1) throw new Exception("Exactly one native physical device required");
            IntPtr physicalList = Allocate(8, owned);
            Check(enumerate(instance, ref physicalCount, physicalList), "physical devices");
            if (physicalCount != 1) throw new Exception("Physical inventory changed");
            IntPtr physical = Marshal.ReadIntPtr(physicalList);
            IntPtr floats = Allocate(88, owned); Marshal.WriteInt32(floats, 1000197000);
            IntPtr driver = Allocate(536, owned); Marshal.WriteInt32(driver, 1000196000); Marshal.WriteIntPtr(driver, 8, floats);
            IntPtr props = Allocate(2048, owned); Marshal.WriteInt32(props, 1000059001); Marshal.WriteIntPtr(props, 8, driver);
            properties(physical, props); IntPtr basic = IntPtr.Add(props, 16);
            string deviceName = Marshal.PtrToStringAnsi(IntPtr.Add(basic, 20));
            string driverName = Marshal.PtrToStringAnsi(IntPtr.Add(driver, 20));
            string driverInfo = Marshal.PtrToStringAnsi(IntPtr.Add(driver, 276));
            if (deviceName != "AMD Radeon 780M Graphics" || driverInfo != "26.8.1 (LLPC)" || U32(basic, 8) != 4098 ||
                U32(basic, 12) != 6400 || U32(driver, 16) != 1 ||
                U32(basic, 0) != 4211037 || U32(floats, 40) != 1 || U32(floats, 64) != 1 ||
                (wide && U32(floats, 68) != 1))
                throw new Exception("Native Radeon identity/retained float controls mismatch");
            bool exactAmdLoaded = false;
            foreach (System.Diagnostics.ProcessModule item in System.Diagnostics.Process.GetCurrentProcess().Modules) {
                if (!Path.GetFileName(item.FileName).Equals("amdvlk64.dll", StringComparison.OrdinalIgnoreCase)) continue;
                if (!item.FileName.Equals(@"C:\WINDOWS\System32\DriverStore\FileRepository\amdvlk.inf_amd64_914ba89eaaafdf60\amdvlk64.dll", StringComparison.OrdinalIgnoreCase) ||
                    Hash(File.ReadAllBytes(item.FileName)) != "5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650")
                    throw new Exception("Loaded AMD ICD path/hash differs from accepted native runtime");
                exactAmdLoaded = true;
            }
            if (!exactAmdLoaded) throw new Exception("Expected native AMD ICD absent");
            if(U32(floats,28)!=1)throw new Exception("SignedZeroInfNanPreserve32 required for exact FMA residual");
            var enumerateExtensions=Bind<DeviceExtensions>(instanceProc(instance,"vkEnumerateDeviceExtensionProperties"),"vkEnumerateDeviceExtensionProperties");
            uint extensionCount=0;
            Check(enumerateExtensions(physical,IntPtr.Zero,ref extensionCount,IntPtr.Zero),"extension count");
            if(extensionCount==0 || extensionCount>1024)throw new Exception("Bounded extension census required");
            uint extensionCapacity=extensionCount;
            IntPtr extensionData=Allocate(checked((int)extensionCapacity*260),owned);
            Check(enumerateExtensions(physical,IntPtr.Zero,ref extensionCount,extensionData),"complete extension census");
            if(extensionCount!=extensionCapacity)throw new Exception("Extension census changed");
            bool fmaExtension=false;
            for(uint i=0;i<extensionCount;++i){
                IntPtr item=IntPtr.Add(extensionData,checked((int)i*260));
                int length=0;while(length<256 && Marshal.ReadByte(item,length)!=0)++length;
                if(length==0 || length==256)throw new Exception("Invalid bounded extension name");
                if(Marshal.PtrToStringAnsi(item,length)=="VK_KHR_shader_fma")fmaExtension=true;
            }
            if(!fmaExtension)throw new Exception("VK_KHR_shader_fma not advertised");
            IntPtr fmaData=Allocate(32,owned);Marshal.WriteInt32(fmaData,1000579000);
            IntPtr features2Data=Allocate(240,owned);Marshal.WriteInt32(features2Data,1000059000);Marshal.WriteIntPtr(features2Data,8,fmaData);
            var features2=Bind<Properties>(instanceProc(instance,"vkGetPhysicalDeviceFeatures2"),"vkGetPhysicalDeviceFeatures2");
            features2(physical,features2Data);
            if(U32(fmaData,20)!=1)throw new Exception("shaderFmaFloat32 not supported");
            Marshal.WriteInt32(fmaData,16,0);Marshal.WriteInt32(fmaData,20,1);Marshal.WriteInt32(fmaData,24,0);
            byte[] extensionName=Encoding.ASCII.GetBytes("VK_KHR_shader_fma\0");
            IntPtr extensionNameData=Allocate(extensionName.Length,owned);Marshal.Copy(extensionName,0,extensionNameData,extensionName.Length);
            IntPtr enabledExtensions=Allocate(8,owned);Marshal.WriteIntPtr(enabledExtensions,extensionNameData);
            checkpoint("fma32-and-signedzero32-supported");
            var getFeatures = Bind<Properties>(instanceProc(instance, "vkGetPhysicalDeviceFeatures"), "vkGetPhysicalDeviceFeatures");
            IntPtr supported = Allocate(220, owned); getFeatures(physical, supported);
            if (U32(supported, 156) != 1) throw new Exception("Native binary64 rung unavailable");
            IntPtr enabled = Allocate(220, owned);
            Marshal.WriteInt32(enabled, 156, 1);


            uint familyCount = 0; families(physical, ref familyCount, IntPtr.Zero);
            if (familyCount == 0 || familyCount > 64) throw new Exception("Invalid queue-family count");
            uint familyCapacity = familyCount; IntPtr familyList = Allocate(checked((int)familyCount * 24), owned);
            families(physical, ref familyCount, familyList);
            if (familyCount == 0 || familyCount > familyCapacity) throw new Exception("Queue inventory grew/truncated");
            uint family = familyCount;
            for (uint i = 0; i < familyCount; ++i)
                if ((U32(familyList, checked((int)i * 24)) & 2) != 0 && U32(familyList, checked((int)i * 24 + 4)) > 0) { family = i; break; }
            if (family == familyCount) throw new Exception("No compute queue");
            IntPtr priority = Allocate(4, owned); Marshal.WriteInt32(priority, 0x3f800000);
            var queueInfo = new QueueInfo { sType = 2, family = family, count = 1, priorities = priority };
            // Match production logical-device Float64/FMA32; modules remain native-default.
            var deviceInfo = new DeviceInfo { sType = 3, pNext = fmaData, queueCount = 1, queues = Structure(queueInfo, owned), features = enabled, extensionCount = 1, extensions = enabledExtensions };
            IntPtr device; Check(createDevice(physical, ref deviceInfo, IntPtr.Zero, out device), "vkCreateDevice");
            var destroyDevice = Bind<DestroyDevice>(deviceProc(device, "vkDestroyDevice"), "vkDestroyDevice");
            cleanup.Add(delegate { destroyDevice(device, IntPtr.Zero); });
            var getQueue = Bind<GetQueue>(deviceProc(device, "vkGetDeviceQueue"), "vkGetDeviceQueue");
            var createBuffer = Bind<CreateObject>(deviceProc(device, "vkCreateBuffer"), "vkCreateBuffer");
            var destroyBuffer = Bind<DestroyObject>(deviceProc(device, "vkDestroyBuffer"), "vkDestroyBuffer");
            var requirements = Bind<BufferRequirements>(deviceProc(device, "vkGetBufferMemoryRequirements"), "vkGetBufferMemoryRequirements");
            var allocateMemory = Bind<CreateObject>(deviceProc(device, "vkAllocateMemory"), "vkAllocateMemory");
            var freeMemory = Bind<DestroyObject>(deviceProc(device, "vkFreeMemory"), "vkFreeMemory");
            var bindMemory = Bind<BindMemory>(deviceProc(device, "vkBindBufferMemory"), "vkBindBufferMemory");
            var mapMemory = Bind<MapMemory>(deviceProc(device, "vkMapMemory"), "vkMapMemory");
            var unmapMemory = Bind<UnmapMemory>(deviceProc(device, "vkUnmapMemory"), "vkUnmapMemory");
            var createShader = Bind<CreateObject>(deviceProc(device, "vkCreateShaderModule"), "vkCreateShaderModule");
            var destroyShader = Bind<DestroyObject>(deviceProc(device, "vkDestroyShaderModule"), "vkDestroyShaderModule");
            var createSetLayout = Bind<CreateObject>(deviceProc(device, "vkCreateDescriptorSetLayout"), "vkCreateDescriptorSetLayout");
            var destroySetLayout = Bind<DestroyObject>(deviceProc(device, "vkDestroyDescriptorSetLayout"), "vkDestroyDescriptorSetLayout");
            var createPipelineLayout = Bind<CreateObject>(deviceProc(device, "vkCreatePipelineLayout"), "vkCreatePipelineLayout");
            var destroyPipelineLayout = Bind<DestroyObject>(deviceProc(device, "vkDestroyPipelineLayout"), "vkDestroyPipelineLayout");
            var createPipelines = Bind<CreatePipelines>(deviceProc(device, "vkCreateComputePipelines"), "vkCreateComputePipelines");
            var destroyPipeline = Bind<DestroyObject>(deviceProc(device, "vkDestroyPipeline"), "vkDestroyPipeline");
            var createPool = Bind<CreateObject>(deviceProc(device, "vkCreateDescriptorPool"), "vkCreateDescriptorPool");
            var destroyPool = Bind<DestroyObject>(deviceProc(device, "vkDestroyDescriptorPool"), "vkDestroyDescriptorPool");
            var allocateSets = Bind<AllocateSets>(deviceProc(device, "vkAllocateDescriptorSets"), "vkAllocateDescriptorSets");
            var updateSets = Bind<UpdateSets>(deviceProc(device, "vkUpdateDescriptorSets"), "vkUpdateDescriptorSets");
            var createCommandPool = Bind<CreateObject>(deviceProc(device, "vkCreateCommandPool"), "vkCreateCommandPool");
            var destroyCommandPool = Bind<DestroyObject>(deviceProc(device, "vkDestroyCommandPool"), "vkDestroyCommandPool");
            var allocateCommands = Bind<AllocateCommands>(deviceProc(device, "vkAllocateCommandBuffers"), "vkAllocateCommandBuffers");
            var begin = Bind<BeginCommand>(deviceProc(device, "vkBeginCommandBuffer"), "vkBeginCommandBuffer");
            var end = Bind<EndCommand>(deviceProc(device, "vkEndCommandBuffer"), "vkEndCommandBuffer");
            var barrier = Bind<Barrier>(deviceProc(device, "vkCmdPipelineBarrier"), "vkCmdPipelineBarrier");
            var bindPipeline = Bind<BindPipeline>(deviceProc(device, "vkCmdBindPipeline"), "vkCmdBindPipeline");
            var bindSets = Bind<BindSets>(deviceProc(device, "vkCmdBindDescriptorSets"), "vkCmdBindDescriptorSets");
            var dispatch = Bind<Dispatch>(deviceProc(device, "vkCmdDispatch"), "vkCmdDispatch");
            var createFence = Bind<CreateObject>(deviceProc(device, "vkCreateFence"), "vkCreateFence");
            var destroyFence = Bind<DestroyObject>(deviceProc(device, "vkDestroyFence"), "vkDestroyFence");
            var submit = Bind<Submit>(deviceProc(device, "vkQueueSubmit"), "vkQueueSubmit");
            var wait = Bind<WaitFences>(deviceProc(device, "vkWaitForFences"), "vkWaitForFences");
            var resetCommand = Bind<ResetCommand>(deviceProc(device, "vkResetCommandBuffer"), "vkResetCommandBuffer");
            var resetFences = Bind<ResetFences>(deviceProc(device, "vkResetFences"), "vkResetFences");
            var queueIdle = Bind<QueueIdle>(deviceProc(device, "vkQueueWaitIdle"), "vkQueueWaitIdle");
            IntPtr queue; getQueue(device, family, 0, out queue);
            if (queue == IntPtr.Zero) throw new Exception("Null compute queue");

            IntPtr memoryPropertiesData = Allocate(520, owned); memoryProperties(physical, memoryPropertiesData);
            uint typeCount = U32(memoryPropertiesData, 0), heapCount = U32(memoryPropertiesData, 260);
            if (typeCount == 0 || typeCount > 32 || heapCount == 0 || heapCount > 16) throw new Exception("Invalid memory counts");
            var buffers = new ulong[2]; var memories = new ulong[2];
            byte[] canary = new byte[expected.Length];
            for (int i = 0; i < canary.Length; ++i) canary[i] = 255;
            byte[][] bytesToWrite = { input, canary }; ulong allocationBytes = 0;
            var allocations = new List<object>();
            for (int row = 0; row < 2; ++row) {
                var info = new BufferInfo { sType = 12, size = (ulong)bytesToWrite[row].Length, usage = 0x20 };
                ulong buffer; Check(createBuffer(device, Structure(info, owned), IntPtr.Zero, out buffer), "vkCreateBuffer");
                buffers[row] = buffer;
                // Register immediately, including requirement/allocation failures.
                // The closure destroys the buffer before freeing bound memory.
                ulong memory = 0;
                cleanup.Add(delegate {
                    destroyBuffer(device, buffer, IntPtr.Zero);
                    if (memory != 0) freeMemory(device, memory, IntPtr.Zero);
                });
                Requirements requirement; requirements(device, buffer, out requirement);
                if (requirement.size < info.size || requirement.alignment == 0 ||
                    requirement.size > AllocationLimit || allocationBytes > AllocationLimit - requirement.size)
                    throw new Exception("Actual Vulkan allocation requirement exceeds diagnostic limit");
                uint selected = typeCount;
                for (uint i = 0; i < typeCount; ++i) {
                    uint flags = U32(memoryPropertiesData, 4 + checked((int)i * 8));
                    uint heap = U32(memoryPropertiesData, 8 + checked((int)i * 8));
                    if (heap >= heapCount) throw new Exception("Invalid memory heap index");
                    if ((requirement.memoryTypeBits & (1u << (int)i)) != 0 && (flags & 6) == 6 &&
                        U64(memoryPropertiesData, 264 + checked((int)heap * 16)) >= requirement.size) { selected = i; break; }
                }
                if (selected == typeCount) throw new Exception("No compatible coherent host-visible allocation");
                var allocation = new AllocationInfo { sType = 5, size = requirement.size, memoryType = selected };
                Check(allocateMemory(device, Structure(allocation, owned), IntPtr.Zero, out memory), "vkAllocateMemory");
                Check(bindMemory(device, buffer, memory, 0), "vkBindBufferMemory");
                memories[row] = memory; allocationBytes += requirement.size;
                allocations.Add(Record("requested_bytes", info.size, "requirement_bytes", requirement.size,
                    "alignment", requirement.alignment, "memory_type", selected,
                    "memory_flags", U32(memoryPropertiesData, 4 + checked((int)selected * 8))));
                IntPtr mapped; Check(mapMemory(device, memory, 0, info.size, 0, out mapped), "vkMapMemory write");
                try { Marshal.Copy(bytesToWrite[row], 0, mapped, bytesToWrite[row].Length); }
                finally { unmapMemory(device, memory); }
            }

            IntPtr code = Allocate(module.Length, owned); Marshal.Copy(module, 0, code, module.Length);
            var shaderInfo = new ShaderInfo { sType = 16, codeSize = new UIntPtr((uint)module.Length), code = code };
            ulong shader; Check(createShader(device, Structure(shaderInfo, owned), IntPtr.Zero, out shader), "vkCreateShaderModule");
            cleanup.Add(delegate { destroyShader(device, shader, IntPtr.Zero); });
            bool verifiedIcd = false;
            foreach (System.Diagnostics.ProcessModule item in System.Diagnostics.Process.GetCurrentProcess().Modules)
                if (Path.GetFileName(item.FileName).Equals("amdvlk64.dll", StringComparison.OrdinalIgnoreCase)) {
                    if (Hash(File.ReadAllBytes(item.FileName)) != "5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650")
                        throw new Exception("Actual loaded AMD ICD changed");
                    verifiedIcd = true;
                }
            if (!verifiedIcd) throw new Exception("Actual native AMD ICD not observed");
            checkpoint("native-device-verified=" + deviceName + ";driver=" + driverInfo + ";denorm32=" + U32(floats, 40) + ";RTE32=" + U32(floats, 64));
            File.WriteAllText(actualPath + ".providers.tsv", "");
            foreach (System.Diagnostics.ProcessModule item in System.Diagnostics.Process.GetCurrentProcess().Modules) {
                string name = Path.GetFileName(item.FileName);
                if (name.Equals("vulkan-1.dll", StringComparison.OrdinalIgnoreCase) || name.Equals("amdvlk64.dll", StringComparison.OrdinalIgnoreCase))
                    File.AppendAllText(actualPath + ".providers.tsv", item.FileName + "\t" + new FileInfo(item.FileName).Length + "\t" + Hash(File.ReadAllBytes(item.FileName)) + Environment.NewLine);
            }
            var bindings = new Binding[] { new Binding { binding = 0, type = 7, count = 1, stages = 0x20 },
                                          new Binding { binding = 1, type = 7, count = 1, stages = 0x20 } };
            var setInfo = new SetLayoutInfo { sType = 32, count = 2, bindings = Structures(bindings, owned) };
            ulong setLayout; Check(createSetLayout(device, Structure(setInfo, owned), IntPtr.Zero, out setLayout), "vkCreateDescriptorSetLayout");
            cleanup.Add(delegate { destroySetLayout(device, setLayout, IntPtr.Zero); });
            IntPtr layouts = Allocate(8, owned); Marshal.WriteInt64(layouts, unchecked((long)setLayout));
            var layoutInfo = new PipelineLayoutInfo { sType = 30, setCount = 1, sets = layouts };
            ulong pipelineLayout; Check(createPipelineLayout(device, Structure(layoutInfo, owned), IntPtr.Zero, out pipelineLayout), "vkCreatePipelineLayout");
            cleanup.Add(delegate { destroyPipelineLayout(device, pipelineLayout, IntPtr.Zero); });
            IntPtr entry = Allocate(5, owned); Marshal.Copy(Encoding.ASCII.GetBytes("main\0"), 0, entry, 5);
            var pipelineInfo = new PipelineInfo { sType = 29, layout = pipelineLayout, baseIndex = -1,
                stage = new ShaderStage { sType = 18, stage = 0x20, module = shader, name = entry } };
            IntPtr pipelinePointer = Allocate(8, owned);
            checkpoint("original-pipeline-begin");
            int pipelineResult = createPipelines(device, 0, 1, Structure(pipelineInfo, owned), IntPtr.Zero, pipelinePointer);
            ulong pipeline = U64(pipelinePointer, 0);
            if (pipeline != 0) cleanup.Add(delegate { destroyPipeline(device, pipeline, IntPtr.Zero); });
            Check(pipelineResult, "vkCreateComputePipelines");
            checkpoint("original-pipeline-complete");
            if (pipeline == 0) throw new Exception("Null compute pipeline");
            IntPtr alternativeCode = Allocate(alternative.Length, owned);
            Marshal.Copy(alternative, 0, alternativeCode, alternative.Length);
            var alternativeShaderInfo = new ShaderInfo { sType = 16, codeSize = new UIntPtr((uint)alternative.Length), code = alternativeCode };
            ulong alternativeShader;
            Check(createShader(device, Structure(alternativeShaderInfo, owned), IntPtr.Zero, out alternativeShader), "vkCreateShaderModule alternative");
            cleanup.Add(delegate { destroyShader(device, alternativeShader, IntPtr.Zero); });
            pipelineInfo.stage.module = alternativeShader;
            IntPtr alternativePointer = Allocate(8, owned);
            checkpoint("candidate-pipeline-begin");
            int alternativeResult = createPipelines(device, 0, 1, Structure(pipelineInfo, owned), IntPtr.Zero, alternativePointer);
            ulong alternativePipeline = U64(alternativePointer, 0);
            if (alternativePipeline != 0) cleanup.Add(delegate { destroyPipeline(device, alternativePipeline, IntPtr.Zero); });
            Check(alternativeResult, "vkCreateComputePipelines alternative");
            checkpoint("candidate-pipeline-complete");
            if (alternativePipeline == 0) throw new Exception("Null alternative compute pipeline");
            var poolInfo = new PoolInfo { sType = 33, maxSets = 1, count = 1,
                sizes = Structure(new PoolSize { type = 7, count = 2 }, owned) };
            ulong pool; Check(createPool(device, Structure(poolInfo, owned), IntPtr.Zero, out pool), "vkCreateDescriptorPool");
            cleanup.Add(delegate { destroyPool(device, pool, IntPtr.Zero); });
            var setAllocation = new SetAllocationInfo { sType = 34, pool = pool, count = 1, layouts = layouts };
            IntPtr set = Allocate(8, owned); Check(allocateSets(device, Structure(setAllocation, owned), set), "vkAllocateDescriptorSets");
            var descriptorBuffers = new DescriptorBuffer[] {
                new DescriptorBuffer { buffer = buffers[0], range = (ulong)input.Length },
                new DescriptorBuffer { buffer = buffers[1], range = (ulong)expected.Length } };
            IntPtr descriptorData = Structures(descriptorBuffers, owned);
            var writes = new WriteSet[] {
                new WriteSet { sType = 35, set = U64(set, 0), binding = 0, count = 1, type = 7, buffers = descriptorData },
                new WriteSet { sType = 35, set = U64(set, 0), binding = 1, count = 1, type = 7, buffers = IntPtr.Add(descriptorData, 24) } };
            updateSets(device, 2, Structures(writes, owned), 0, IntPtr.Zero);
            var commandPoolInfo = new CommandPoolInfo { sType = 39, flags = 2, family = family };
            ulong commandPool; Check(createCommandPool(device, Structure(commandPoolInfo, owned), IntPtr.Zero, out commandPool), "vkCreateCommandPool");
            cleanup.Add(delegate { destroyCommandPool(device, commandPool, IntPtr.Zero); });
            var commandAllocation = new CommandAllocationInfo { sType = 40, pool = commandPool, count = 1 };
            IntPtr command; Check(allocateCommands(device, Structure(commandAllocation, owned), out command), "vkAllocateCommandBuffers");
            ulong fence; Check(createFence(device, Structure(new FenceInfo { sType = 8 }, owned), IntPtr.Zero, out fence), "vkCreateFence");
            cleanup.Add(delegate { destroyFence(device, fence, IntPtr.Zero); });
            IntPtr commands = Allocate(8, owned); Marshal.WriteIntPtr(commands, command);
            IntPtr fences = Allocate(8, owned); Marshal.WriteInt64(fences, unchecked((long)fence));
            var submitInfo = new SubmitInfo { sType = 4, commandCount = 1, commands = commands };
            var samples = new List<object>(); byte[] actual = new byte[expected.Length];
            byte[] inputShadow = (byte[])input.Clone();
            if (sequencePaths.Length == 0 || sequencePaths.Length > 32 ||
                sequenceBindings.Length != sequencePaths.Length || activeCounts.Length != sequencePaths.Length ||
                uploadBytes.Length != sequencePaths.Length || sequenceNames.Length != sequencePaths.Length ||
                outputWords != checked((int)Rows * 4699)) throw new Exception("Sequence shape mismatch");
            uint dispatches = 0;
            for (int sequence = 0; sequence < sequencePaths.Length; ++sequence) {
                uint active = activeCounts[sequence];
                if (active == 0 || active > Rows ||
                    !System.Text.RegularExpressions.Regex.IsMatch(sequenceNames[sequence], "^[a-z0-9-]+$"))
                    throw new Exception("Active prefix/name invalid");
                byte[] update = ReadBounded(sequencePaths[sequence]);
                if (Hash(update) != sequenceBindings[sequence] || update.Length != uploadBytes[sequence] ||
                    Word(update, 0) != Rows ||
                    (update.Length != checked(4 + (int)active * 920) && update.Length != input.Length))
                    throw new Exception("Exact uploaded prefix/table binding mismatch");
                Array.Copy(update, 0, inputShadow, 0, update.Length);
                IntPtr uploaded; Check(mapMemory(device, memories[0], 0, (ulong)update.Length, 0, out uploaded), "vkMapMemory prefix");
                try { Marshal.Copy(update, 0, uploaded, update.Length); }
                finally { unmapMemory(device, memories[0]); }
                byte[] baseline = null;
                for (int variant = 0; variant < 2; ++variant) {
                    bool candidate = variant == 1;
                    IntPtr clear; Check(mapMemory(device, memories[1], 0, (ulong)actual.Length, 0, out clear), "vkMapMemory canary");
                    try { Marshal.Copy(canary, 0, clear, canary.Length); }
                    finally { unmapMemory(device, memories[1]); }
                    Check(resetCommand(command, 0), "vkResetCommandBuffer");
                    Check(resetFences(device, 1, fences), "vkResetFences");
                    var beginInfo = new CommandBeginInfo { sType = 42, flags = 1 };
                    Check(begin(command, ref beginInfo), "vkBeginCommandBuffer");
                    var before = new MemoryBarrier { sType = 46, source = 0x40 | 0x4000, destination = 0x20 | 0x40 };
                    barrier(command, 0x800 | 0x4000, 0x800, 0, 1, ref before, 0, IntPtr.Zero, 0, IntPtr.Zero);
                    bindPipeline(command, 1, candidate ? alternativePipeline : pipeline);
                    bindSets(command, 1, pipelineLayout, 0, 1, set, 0, IntPtr.Zero);
                    uint groupsX = candidate ? (active + 1) / 2 : active;
                    uint groupsY = candidate && (active & 1) != 0 ? 2u : 1u;
                    dispatch(command, groupsX, groupsY, 1);
                    var after = new MemoryBarrier { sType = 46, source = 0x40, destination = 0x20 | 0x2000 };
                    barrier(command, 0x800, 0x800 | 0x4000, 0, 1, ref after, 0, IntPtr.Zero, 0, IntPtr.Zero);
                    Check(end(command), "vkEndCommandBuffer");
                    submitted = false; completed = false;
                    Check(submit(queue, 1, ref submitInfo, fence), "vkQueueSubmit"); submitted = true;
                    Check(wait(device, 1, fences, 1, FenceTimeoutNanoseconds), "bounded vkWaitForFences");
                    completed = true; ++dispatches;
                    IntPtr read; Check(mapMemory(device, memories[1], 0, (ulong)actual.Length, 0, out read), "vkMapMemory read");
                    try { Marshal.Copy(read, actual, 0, actual.Length); }
                    finally { unmapMemory(device, memories[1]); }
                    byte[] inputAfter = new byte[inputShadow.Length];
                    Check(mapMemory(device, memories[0], 0, (ulong)inputAfter.Length, 0, out read), "vkMapMemory input seal");
                    try { Marshal.Copy(read, inputAfter, 0, inputAfter.Length); }
                    finally { unmapMemory(device, memories[0]); }
                    if (Hash(inputAfter) != Hash(inputShadow)) throw new Exception("Device changed input/table words");
                    string path = actualPath + "." + sequenceNames[sequence] + (candidate ? ".candidate.bin" : ".baseline.bin");
                    File.WriteAllBytes(path, actual);
                    if (baseline == null) baseline = (byte[])actual.Clone();
                    uint differences = 0, inactiveChanges = 0;
                    var examples = new List<object>();
                    for (int i = 0; i < outputWords; ++i) {
                        uint a = Word(actual, i * 4), e = Word(baseline, i * 4);
                        if (a != e) { ++differences; if (examples.Count < 8) examples.Add(Record("word", i, "expected", e, "actual", a)); }
                        if ((uint)i >= active * 4699 && a != 0xffffffffu) ++inactiveChanges;
                    }
                    samples.Add(Record("sequence", sequenceNames[sequence], "variant", candidate ? "candidate" : "baseline",
                        "active_rows", active, "groups_x", groupsX, "groups_y", groupsY,
                        "uploaded_bytes", update.Length, "input_shadow_sha256", Hash(inputShadow),
                        "output_file", Path.GetFileName(path), "output_sha256", Hash(actual),
                        "whole_capacity_word_mismatches", differences, "inactive_word_changes", inactiveChanges));
                    checkpoint("sequence=" + sequenceNames[sequence] + ";variant=" + variant + ";differences=" + differences + ";inactive=" + inactiveChanges);
                    File.AppendAllText(actualPath + ".observations.tsv", String.Join("\t", new string[] {
                        sequenceNames[sequence], variant.ToString(), active.ToString(), groupsX.ToString(), groupsY.ToString(),
                        update.Length.ToString(), Hash(inputShadow), Hash(actual), differences.ToString(), inactiveChanges.ToString()
                    }) + Environment.NewLine);
                    if (differences != 0 || inactiveChanges != 0) throw new Exception("Paired prefix changed complete output or inactive canary; examples=" + examples.Count);
                }
            }
            var loadedModules = new List<object>();
            foreach (System.Diagnostics.ProcessModule item in System.Diagnostics.Process.GetCurrentProcess().Modules) {
                string name = Path.GetFileName(item.FileName);
                if (name.Equals("vulkan-1.dll", StringComparison.OrdinalIgnoreCase) || name.Equals("amdvlk64.dll", StringComparison.OrdinalIgnoreCase))
                    loadedModules.Add(Record("path", item.FileName, "bytes", new FileInfo(item.FileName).Length, "sha256", Hash(File.ReadAllBytes(item.FileName))));
            }
            return Record("completed", true, "scope", "Isolated paired-row raw/reference gate; no timing or full science/release acceptance",
                "device", deviceName, "driver", driverInfo, "shader_fma32_enabled", true, "shader_float64_enabled", true,
                "loaded_modules", loadedModules, "allocations", allocations, "actual_resident_bytes", allocationBytes,
                "module_sha256", moduleHash, "alternative_module_sha256", Hash(alternative),
                "capacity", Rows, "completed_dispatches", dispatches, "fence_timeout_ns", FenceTimeoutNanoseconds,
                "observed_words_per_dispatch", outputWords, "samples", samples);
        } finally {
            if (submitted && !completed) {
                // In-flight resources must not be destroyed or the DLL unloaded.
                // The isolated child exits immediately; its native supervisor
                // bounds termination even if the driver does not return.
                Console.Error.WriteLine("Fence not complete: native resources retained until isolated process termination.");
            } else {
                for (int i = cleanup.Count - 1; i >= 0; --i) cleanup[i]();
                foreach (IntPtr pointer in owned) Marshal.FreeHGlobal(pointer);
                FreeLibrary(loader);
            }
        }
    }
}
