// Baseline-only ordinary native pipeline generated-code observer.
// Adds only the advertised AMD query extension; no dispatch or qualification claim.
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;

public static class NativeTransportFmaProbe {
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


    [StructLayout(LayoutKind.Sequential)] struct CacheInfo {
        public uint sType; public IntPtr pNext; public uint flags; public UIntPtr initialBytes; public IntPtr initialData;
    }
    [UnmanagedFunctionPointer(CallingConvention.Winapi)] delegate int ShaderQuery(IntPtr device, ulong pipeline, uint stage, uint type, ref UIntPtr size, IntPtr data);
    static object QueryOne(ShaderQuery query, IntPtr device, ulong pipeline, uint type, string path) {
        UIntPtr size = UIntPtr.Zero;
        int first = query(device, pipeline, 0x20, type, ref size, IntPtr.Zero);
        if (first == -8) return Record("type", type, "available", false, "size_query_result", first, "reason", "VK_ERROR_FEATURE_NOT_PRESENT");
        Check(first, "AMD shader-info size query");
        ulong capacity = size.ToUInt64();
        if (capacity == 0 || capacity > 8 * 1024 * 1024) throw new Exception("Invalid bounded shader-info size");
        IntPtr data = Marshal.AllocHGlobal(checked((int)capacity));
        try {
            for (int i = 0; i < (int)capacity; ++i) Marshal.WriteByte(data, i, 0);
            int complete = query(device, pipeline, 0x20, type, ref size, data);
            Check(complete, "AMD shader-info complete query; no truncated acceptance");
            if (size.ToUInt64() == 0 || size.ToUInt64() > capacity) throw new Exception("Shader-info query exceeded capacity");
            byte[] bytes = new byte[checked((int)size.ToUInt64())]; Marshal.Copy(data, bytes, 0, bytes.Length);
            if (type == 2 && bytes[bytes.Length - 1] != 0) throw new Exception("Disassembly lacks required terminal zero");
            File.WriteAllBytes(path, bytes);
            object statistics = null;
            if (type == 0) {
                if (bytes.Length != 72) throw new Exception("Unexpected Win64 AMD statistics size");
                statistics = Record("shader_stage_mask", U32(data, 0), "used_vgprs", U32(data, 8), "used_sgprs", U32(data, 12),
                    "lds_bytes_per_workgroup_limit", U32(data, 16), "lds_usage_bytes", U64(data, 24), "scratch_usage_bytes", U64(data, 32),
                    "physical_vgprs", U32(data, 40), "physical_sgprs", U32(data, 44), "available_vgprs", U32(data, 48), "available_sgprs", U32(data, 52),
                    "workgroup_x", U32(data, 56), "workgroup_y", U32(data, 60), "workgroup_z", U32(data, 64));
                if (U32(data, 0) != 0x20 || U32(data, 56) != 64 || U32(data, 60) != 1 || U32(data, 64) != 1) throw new Exception("Unexpected compute statistics identity");
            }
            return Record("type", type, "available", true, "size_query_result", first, "complete_query_result", complete,
                "capacity_bytes", capacity, "returned_bytes", bytes.Length, "sha256", Hash(bytes), "file", Path.GetFileName(path), "statistics", statistics);
        } finally { Marshal.FreeHGlobal(data); }
    }
    public static object Observe(string loaderPath, string transportPath, string endpointPath, string outputDirectory) {
        VerifyAbi(); Abi<CacheInfo>(40);
        if (Hash(ReadBounded(loaderPath)) != DllHash) throw new Exception("Pinned loader hash mismatch");
        IntPtr loader = LoadLibraryEx(loaderPath, IntPtr.Zero, 8);
        if (loader == IntPtr.Zero) throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error());
        var owned = new List<IntPtr>(); var cleanup = new List<Action>();
        try {
            var instanceProc = Bind<ProcAddress>(GetProcAddress(loader, "vkGetInstanceProcAddr"), "vkGetInstanceProcAddr");
            var createInstance = Bind<CreateInstance>(instanceProc(IntPtr.Zero, "vkCreateInstance"), "vkCreateInstance");
            var app = new ApplicationInfo { sType = 0, apiVersion = (1u << 22) | (3u << 12) };
            var info = new InstanceInfo { sType = 1, application = Structure(app, owned) };
            IntPtr instance; Check(createInstance(ref info, IntPtr.Zero, out instance), "vkCreateInstance");
            var destroy = Bind<DestroyInstance>(instanceProc(instance, "vkDestroyInstance"), "vkDestroyInstance");
            cleanup.Add(delegate { destroy(instance, IntPtr.Zero); });
            var enumerate = Bind<EnumerateDevices>(instanceProc(instance, "vkEnumeratePhysicalDevices"), "vkEnumeratePhysicalDevices");
            uint count = 0; Check(enumerate(instance, ref count, IntPtr.Zero), "physical count");
            if (count != 1) throw new Exception("Exactly one native physical device required");
            IntPtr data = Allocate(8, owned); Check(enumerate(instance, ref count, data), "complete physical census");
            if (count != 1) throw new Exception("Physical census changed");
            IntPtr physical = Marshal.ReadIntPtr(data);
            var properties = Bind<Properties>(instanceProc(instance, "vkGetPhysicalDeviceProperties2"), "vkGetPhysicalDeviceProperties2");
            IntPtr floats = Allocate(88, owned); Marshal.WriteInt32(floats, 1000197000);
            IntPtr driver = Allocate(536, owned); Marshal.WriteInt32(driver, 1000196000); Marshal.WriteIntPtr(driver, 8, floats);
            IntPtr props = Allocate(2048, owned); Marshal.WriteInt32(props, 1000059001); Marshal.WriteIntPtr(props, 8, driver);
            properties(physical, props); IntPtr basic = IntPtr.Add(props, 16);
            string deviceName = Marshal.PtrToStringAnsi(IntPtr.Add(basic, 20));
            string driverInfo = Marshal.PtrToStringAnsi(IntPtr.Add(driver, 276));
            if (deviceName != "AMD Radeon 780M Graphics" || driverInfo != "26.8.1 (LLPC)" || U32(basic, 8) != 4098 || U32(basic, 12) != 6400)
                throw new Exception("Native provider identity mismatch");
            var enumerateExtensions = Bind<DeviceExtensions>(instanceProc(instance, "vkEnumerateDeviceExtensionProperties"), "vkEnumerateDeviceExtensionProperties");
            uint extensionCount = 0; Check(enumerateExtensions(physical, IntPtr.Zero, ref extensionCount, IntPtr.Zero), "extension count");
            if (extensionCount == 0 || extensionCount > 1024) throw new Exception("Unbounded extension census");
            uint capacity = extensionCount; IntPtr extensions = Allocate(checked((int)capacity * 260), owned);
            Check(enumerateExtensions(physical, IntPtr.Zero, ref extensionCount, extensions), "full extension census");
            if (extensionCount != capacity) throw new Exception("Extension census changed");
            var names = new List<object>(); var advertised = new HashSet<string>(); bool amdInfo = false;
            for (uint i = 0; i < extensionCount; ++i) {
                IntPtr item = IntPtr.Add(extensions, checked((int)i * 260));
                int length = 0; while (length < 256 && Marshal.ReadByte(item, length) != 0) ++length;
                if (length == 0 || length == 256) throw new Exception("Invalid extension name");
                string name = Marshal.PtrToStringAnsi(item, length);
                names.Add(Record("name", name, "spec_version", U32(item, 256))); advertised.Add(name);
                if (name == "VK_AMD_shader_info") amdInfo = true;
            }

            if (!amdInfo || !advertised.Contains("VK_KHR_shader_fma")) throw new Exception("Expected query/FMA extensions absent");
            if (advertised.Contains("VK_KHR_portability_subset")) throw new Exception("Unexpected portability device; observer is provider-specific");
            if (U32(floats, 28) != 1 || U32(floats, 40) != 1 || U32(floats, 64) != 1) throw new Exception("Required binary32 controls absent");
            IntPtr fma = Allocate(32, owned); Marshal.WriteInt32(fma, 1000579000);
            IntPtr features2 = Allocate(240, owned); Marshal.WriteInt32(features2, 1000059000); Marshal.WriteIntPtr(features2, 8, fma);
            Bind<Properties>(instanceProc(instance, "vkGetPhysicalDeviceFeatures2"), "vkGetPhysicalDeviceFeatures2")(physical, features2);
            if (U32(fma, 20) != 1 || U32(features2, 16 + 156) != 1) throw new Exception("Required shaderFma32/Float64 features absent");
            Marshal.WriteInt32(fma, 16, 0); Marshal.WriteInt32(fma, 20, 1); Marshal.WriteInt32(fma, 24, 0);
            IntPtr enabled = Allocate(220, owned); Marshal.WriteInt32(enabled, 156, 1);
            string[] enabledNames = new string[] { "VK_KHR_shader_fma", "VK_AMD_shader_info" };
            IntPtr enabledExtensions = Allocate(16, owned);
            for (int i = 0; i < enabledNames.Length; ++i) {
                byte[] name = Encoding.ASCII.GetBytes(enabledNames[i] + "\0");
                IntPtr p = Allocate(name.Length, owned); Marshal.Copy(name, 0, p, name.Length); Marshal.WriteIntPtr(enabledExtensions, i * 8, p);
            }
            var families = Bind<Families>(instanceProc(instance, "vkGetPhysicalDeviceQueueFamilyProperties"), "vkGetPhysicalDeviceQueueFamilyProperties");
            uint familyCount = 0; families(physical, ref familyCount, IntPtr.Zero);
            if (familyCount == 0 || familyCount > 64) throw new Exception("Invalid bounded queue-family count");
            uint familyCapacity = familyCount; IntPtr familyData = Allocate(checked((int)familyCount * 24), owned);
            families(physical, ref familyCount, familyData);
            if (familyCount != familyCapacity) throw new Exception("Queue-family census changed");
            uint family = familyCount;
            for (uint i = 0; i < familyCount; ++i) if ((U32(familyData, (int)i * 24) & 2) != 0 && U32(familyData, (int)i * 24 + 4) > 0) { family = i; break; }
            if (family == familyCount) throw new Exception("No compute queue");
            IntPtr priority = Allocate(4, owned); Marshal.WriteInt32(priority, 0x3f800000);
            var queueInfo = new QueueInfo { sType = 2, family = family, count = 1, priorities = priority };
            var deviceInfo = new DeviceInfo { sType = 3, pNext = fma, queueCount = 1, queues = Structure(queueInfo, owned), features = enabled, extensionCount = 2, extensions = enabledExtensions };
            IntPtr device;
            Check(Bind<CreateDevice>(instanceProc(instance, "vkCreateDevice"), "vkCreateDevice")(physical, ref deviceInfo, IntPtr.Zero, out device), "vkCreateDevice");
            var deviceProc = Bind<ProcAddress>(instanceProc(instance, "vkGetDeviceProcAddr"), "vkGetDeviceProcAddr");
            var destroyDevice = Bind<DestroyDevice>(deviceProc(device, "vkDestroyDevice"), "vkDestroyDevice");
            cleanup.Add(delegate { destroyDevice(device, IntPtr.Zero); });
            var createCache = Bind<CreateObject>(deviceProc(device, "vkCreatePipelineCache"), "vkCreatePipelineCache");
            var destroyCache = Bind<DestroyObject>(deviceProc(device, "vkDestroyPipelineCache"), "vkDestroyPipelineCache");
            ulong cache; Check(createCache(device, Structure(new CacheInfo { sType = 17 }, owned), IntPtr.Zero, out cache), "empty vkCreatePipelineCache");
            cleanup.Add(delegate { destroyCache(device, cache, IntPtr.Zero); });
            var createSet = Bind<CreateObject>(deviceProc(device, "vkCreateDescriptorSetLayout"), "vkCreateDescriptorSetLayout");
            var destroySet = Bind<DestroyObject>(deviceProc(device, "vkDestroyDescriptorSetLayout"), "vkDestroyDescriptorSetLayout");
            var setInfo = new SetLayoutInfo { sType = 32, count = 2, bindings = Structures(new Binding[] {
                new Binding { binding = 0, type = 7, count = 1, stages = 0x20 }, new Binding { binding = 1, type = 7, count = 1, stages = 0x20 } }, owned) };
            ulong set; Check(createSet(device, Structure(setInfo, owned), IntPtr.Zero, out set), "vkCreateDescriptorSetLayout");
            cleanup.Add(delegate { destroySet(device, set, IntPtr.Zero); });
            IntPtr sets = Allocate(8, owned); Marshal.WriteInt64(sets, unchecked((long)set));
            var createLayout = Bind<CreateObject>(deviceProc(device, "vkCreatePipelineLayout"), "vkCreatePipelineLayout");
            var destroyLayout = Bind<DestroyObject>(deviceProc(device, "vkDestroyPipelineLayout"), "vkDestroyPipelineLayout");
            ulong layout; Check(createLayout(device, Structure(new PipelineLayoutInfo { sType = 30, setCount = 1, sets = sets }, owned), IntPtr.Zero, out layout), "vkCreatePipelineLayout");
            cleanup.Add(delegate { destroyLayout(device, layout, IntPtr.Zero); });
            var createShader = Bind<CreateObject>(deviceProc(device, "vkCreateShaderModule"), "vkCreateShaderModule");
            var destroyShader = Bind<DestroyObject>(deviceProc(device, "vkDestroyShaderModule"), "vkDestroyShaderModule");
            var createPipeline = Bind<CreatePipelines>(deviceProc(device, "vkCreateComputePipelines"), "vkCreateComputePipelines");
            var destroyPipeline = Bind<DestroyObject>(deviceProc(device, "vkDestroyPipeline"), "vkDestroyPipeline");
            var query = Bind<ShaderQuery>(deviceProc(device, "vkGetShaderInfoAMD"), "vkGetShaderInfoAMD");
            IntPtr entry = Allocate(5, owned); Marshal.Copy(Encoding.ASCII.GetBytes("main\0"), 0, entry, 5);
            string[] paths = new string[] { transportPath };
            string[] stages = new string[] { "TransportPaired" };
            string[] expectedHashes = new string[] { "721be7625f8eebbf722adf3e79a89050bffe9bf9b7260b8f5120b0a85ccd7f27" };
            var observations = new List<object>();
            for (int i = 0; i < paths.Length; ++i) {
                byte[] bytes = ReadBounded(paths[i]); if (Hash(bytes) != expectedHashes[i]) throw new Exception("Exact module binding mismatch");
                IntPtr code = Allocate(bytes.Length, owned); Marshal.Copy(bytes, 0, code, bytes.Length);
                ulong shader; Check(createShader(device, Structure(new ShaderInfo { sType = 16, codeSize = new UIntPtr((uint)bytes.Length), code = code }, owned), IntPtr.Zero, out shader), "vkCreateShaderModule");
                ulong capturedShader = shader; cleanup.Add(delegate { destroyShader(device, capturedShader, IntPtr.Zero); });
                var pipelineInfo = new PipelineInfo { sType = 29, layout = layout, stage = new ShaderStage { sType = 18, stage = 0x20, module = shader, name = entry } };
                IntPtr result = Allocate(8, owned); var clock = System.Diagnostics.Stopwatch.StartNew();
                Console.WriteLine(stages[i] + " ordinary pipeline begin"); Console.Out.Flush();
                int rc = createPipeline(device, cache, 1, Structure(pipelineInfo, owned), IntPtr.Zero, result); clock.Stop();
                ulong pipeline = U64(result, 0);
                if (pipeline != 0) { ulong capturedPipeline = pipeline; cleanup.Add(delegate { destroyPipeline(device, capturedPipeline, IntPtr.Zero); }); }
                Check(rc, "ordinary vkCreateComputePipelines"); if (pipeline == 0) throw new Exception("Null pipeline");
                var infos = new List<object>();
                for (uint type = 0; type < 3; ++type) infos.Add(QueryOne(query, device, pipeline, type, Path.Combine(outputDirectory, stages[i].ToLowerInvariant() + ".info" + type)));
                observations.Add(Record("stage", stages[i], "module_bytes", bytes.Length, "module_sha256", Hash(bytes), "pipeline_flags", pipelineInfo.flags,
                    "stage_flags", pipelineInfo.stage.flags, "pipeline_creation_result", rc, "pipeline_creation_ms", clock.Elapsed.TotalMilliseconds, "queries", infos));
                Console.WriteLine(stages[i] + " ordinary queries complete"); Console.Out.Flush();
            }
            var modules = new List<object>(); bool amdLoaded = false, loaderLoaded = false;
            foreach (System.Diagnostics.ProcessModule m in System.Diagnostics.Process.GetCurrentProcess().Modules) {
                string name = Path.GetFileName(m.FileName);
                if (!name.Equals("amdvlk64.dll", StringComparison.OrdinalIgnoreCase) && !name.Equals("vulkan-1.dll", StringComparison.OrdinalIgnoreCase)) continue;
                string hash = Hash(File.ReadAllBytes(m.FileName));
                if (name.Equals("amdvlk64.dll", StringComparison.OrdinalIgnoreCase)) {
                    if (!m.FileName.Equals(@"C:\WINDOWS\System32\DriverStore\FileRepository\amdvlk.inf_amd64_914ba89eaaafdf60\amdvlk64.dll", StringComparison.OrdinalIgnoreCase) || hash != "5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650") throw new Exception("Actual AMD ICD differs");
                    amdLoaded = true;
                } else {
                    if (!Path.GetFullPath(m.FileName).Equals(Path.GetFullPath(loaderPath), StringComparison.OrdinalIgnoreCase) || hash != DllHash) throw new Exception("Actual loader differs");
                    loaderLoaded = true;
                }
                modules.Add(Record("name", name, "path", m.FileName, "bytes", new FileInfo(m.FileName).Length, "sha256", hash));
            }
            if (!amdLoaded || !loaderLoaded || modules.Count != 2) throw new Exception("Expected complete actual provider census");

            return Record("completed", true, "instance_api", "1.3", "device", deviceName, "driver", driverInfo, "loaded_modules", modules,
                "shader_float64_enabled", true, "shader_fma32_enabled", true, "shader_fma16_enabled", false, "shader_fma64_enabled", false,
                "enabled_extensions", enabledNames, "queue_family", family, "signed_zero_infnan32", U32(floats, 28), "denorm32", U32(floats, 40), "RTE32", U32(floats, 64),
                "pipeline_cache_flags", 0, "pipeline_cache_initial_bytes", 0, "driver_cold_verified", false, "pipelines_created", 1, "dispatches", 0, "observations", observations,
                "scope", "Isolated unadopted paired-row Transport module; production features/layout/flags, tooling extension added. Static generated-code only; zero dispatches, no numerical, speed, or qualification verdict.");
        } finally {
            for (int i = cleanup.Count - 1; i >= 0; --i) cleanup[i]();
            foreach (IntPtr p in owned) Marshal.FreeHGlobal(p);
            FreeLibrary(loader);
        }
    }
}
