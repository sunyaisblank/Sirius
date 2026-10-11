// Finite native Transport comparison using the existing private Win64 Vulkan harness.
// Reads exact validated modules; this is not complete product qualification.
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

    public static object Census(string loaderPath) {
        VerifyAbi();
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
            var names = new List<object>(); bool amdInfo = false;
            for (uint i = 0; i < extensionCount; ++i) {
                IntPtr item = IntPtr.Add(extensions, checked((int)i * 260));
                int length = 0; while (length < 256 && Marshal.ReadByte(item, length) != 0) ++length;
                if (length == 0 || length == 256) throw new Exception("Invalid extension name");
                string name = Marshal.PtrToStringAnsi(item, length);
                names.Add(Record("name", name, "spec_version", U32(item, 256)));
                if (name == "VK_AMD_shader_info") amdInfo = true;
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
            return Record("completed", true, "instance_api", "1.3", "device", deviceName, "driver", driverInfo,
                "extension_count", extensionCount, "extensions", names, "VK_AMD_shader_info_advertised", amdInfo,
                "shaderFloat64_supported", U32Feature(instanceProc, instance, physical, owned), "loaded_modules", modules,
                "logical_devices_created", 0, "pipelines_created", 0, "dispatches", 0,
                "scope", "Current native extension census only; no generated-code, performance, product or qualification verdict.");
        } finally {
            for (int i = cleanup.Count - 1; i >= 0; --i) cleanup[i]();
            foreach (IntPtr p in owned) Marshal.FreeHGlobal(p);
            FreeLibrary(loader);
        }
    }
    static uint U32Feature(ProcAddress proc, IntPtr instance, IntPtr physical, List<IntPtr> owned) {
        IntPtr data = Allocate(220, owned);
        Bind<Properties>(proc(instance, "vkGetPhysicalDeviceFeatures"), "vkGetPhysicalDeviceFeatures")(physical, data);
        return U32(data, 156);
    }
}
