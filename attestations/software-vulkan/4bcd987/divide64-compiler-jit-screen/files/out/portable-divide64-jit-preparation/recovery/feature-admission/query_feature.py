"""Read-only Vulkan core feature query: no surfaces, device, pipeline or dispatch."""
import ctypes as c
import hashlib
import json
import os
from pathlib import Path
import re
import time

root = Path.cwd()
work = Path(__file__).resolve().parent
assert os.environ['VK_DRIVER_FILES'] == '/usr/share/vulkan/icd.d/lvp_icd.json'
assert not os.environ.get('LD_PRELOAD') and 'VK_ICD_FILENAMES' not in os.environ
paths = ['/usr/include/vulkan/vulkan_core.h', '/usr/share/vulkan/icd.d/lvp_icd.json',
         '/usr/lib/x86_64-linux-gnu/libvulkan.so.1', '/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so',
         '/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1', str(Path(__file__).resolve())]
def pin(name):
    p = Path(name); raw = p.read_bytes()
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(), 'resolved_path': str(p.resolve(strict=True))}
before = {p: pin(p) for p in paths}
header = Path(paths[0]).read_text()
body = re.search(r'typedef struct VkPhysicalDeviceFeatures \{(.*?)\} VkPhysicalDeviceFeatures;', header, re.S)[1]
names = re.findall(r'VkBool32\s+(\w+);', body)
assert len(names) == 55 and names[40] == 'shaderInt64'
class Features(c.Structure):
    _fields_ = [(name, c.c_uint32) for name in names]
class Application(c.Structure):
    _fields_ = [('sType', c.c_uint32), ('pNext', c.c_void_p), ('pApplicationName', c.c_char_p),
                ('applicationVersion', c.c_uint32), ('pEngineName', c.c_char_p),
                ('engineVersion', c.c_uint32), ('apiVersion', c.c_uint32)]
class InstanceInfo(c.Structure):
    _fields_ = [('sType', c.c_uint32), ('pNext', c.c_void_p), ('flags', c.c_uint32),
                ('pApplicationInfo', c.POINTER(Application)), ('enabledLayerCount', c.c_uint32),
                ('ppEnabledLayerNames', c.c_void_p), ('enabledExtensionCount', c.c_uint32),
                ('ppEnabledExtensionNames', c.c_void_p)]
vk = c.CDLL(paths[2])
vk.vkCreateInstance.argtypes = [c.POINTER(InstanceInfo), c.c_void_p, c.POINTER(c.c_void_p)]
vk.vkCreateInstance.restype = c.c_int32
vk.vkEnumeratePhysicalDevices.argtypes = [c.c_void_p, c.POINTER(c.c_uint32), c.POINTER(c.c_void_p)]
vk.vkEnumeratePhysicalDevices.restype = c.c_int32
vk.vkGetPhysicalDeviceFeatures.argtypes = [c.c_void_p, c.POINTER(Features)]
vk.vkGetPhysicalDeviceFeatures.restype = None
vk.vkGetPhysicalDeviceProperties.argtypes = [c.c_void_p, c.c_void_p]
vk.vkGetPhysicalDeviceProperties.restype = None
vk.vkDestroyInstance.argtypes = [c.c_void_p, c.c_void_p]
vk.vkDestroyInstance.restype = None
app = Application(sType=0, pApplicationName=b'Sirius read-only feature query', apiVersion=(1 << 22) | (2 << 12))
info = InstanceInfo(sType=1, pApplicationInfo=c.pointer(app))
instance = c.c_void_p(); started = time.monotonic()
print('vkCreateInstance begin', flush=True)
assert vk.vkCreateInstance(c.byref(info), None, c.byref(instance)) == 0
print('vkCreateInstance complete', flush=True)
try:
    count = c.c_uint32()
    assert vk.vkEnumeratePhysicalDevices(instance, c.byref(count), None) == 0 and count.value == 1
    devices = (c.c_void_p * count.value)()
    assert vk.vkEnumeratePhysicalDevices(instance, c.byref(count), devices) == 0 and count.value == 1
    features = Features(); vk.vkGetPhysicalDeviceFeatures(devices[0], c.byref(features))
    # Core Properties starts with five uint32s, deviceName[256], UUID[16].
    # Overallocate aligned storage for the complete implementation-written core struct.
    properties = (c.c_uint64 * 512)()
    vk.vkGetPhysicalDeviceProperties(devices[0], c.byref(properties))
    raw = bytes(properties); device_name = raw[20:276].split(b'\0', 1)[0].decode()
    assert device_name.startswith('llvmpipe') and c.c_uint32.from_buffer_copy(raw[16:20]).value == 4
    maps = Path('/proc/self/maps').read_text(); (work / 'direct-query-maps.log').write_text(maps)
    for name in paths[2:5]:
        assert str(Path(name).resolve()) in maps, name
    result = {'schema': 'sirius-direct-selected-provider-feature-query-v1', 'physical_devices': count.value,
              'device_name': device_name, 'device_type': 4, 'pipeline_cache_uuid': raw[276:292].hex(),
              'features': {name: int(getattr(features, name)) for name in names},
              'features_struct_bytes': c.sizeof(features), 'header_order_verified': True,
              'no_surface_device_pipeline_or_dispatch_created': True,
              'application_shaderInt64_remains_disabled': True, 'qualification_claimed': False}
finally:
    vk.vkDestroyInstance(instance, None)
result.update(wall_seconds=time.monotonic()-started, before=before, after={p: pin(p) for p in paths},
              maps=pin(work / 'direct-query-maps.log'))
assert result['before'] == result['after']
(work / 'direct-query.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('device_name','pipeline_cache_uuid','wall_seconds','qualification_claimed')},
                 sort_keys=True), flush=True)
print('shaderInt64=' + str(features.shaderInt64), flush=True)
