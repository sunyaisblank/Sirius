#!/usr/bin/env python3
"""Finite public-Vulkan property query; no device creation, pipeline or dispatch."""
import argparse
import ctypes as C
import hashlib
import json
from pathlib import Path
import platform
import sys


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path):
    path = Path(path).resolve(strict=True)
    before = path.stat()
    data = path.read_bytes()
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), "image changed while hashing")
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


class ApplicationInfo(C.Structure):
    _fields_ = [("sType", C.c_uint32), ("pNext", C.c_void_p), ("pApplicationName", C.c_char_p),
                ("applicationVersion", C.c_uint32), ("pEngineName", C.c_char_p),
                ("engineVersion", C.c_uint32), ("apiVersion", C.c_uint32)]


class InstanceInfo(C.Structure):
    _fields_ = [("sType", C.c_uint32), ("pNext", C.c_void_p), ("flags", C.c_uint32),
                ("pApplicationInfo", C.POINTER(ApplicationInfo)), ("enabledLayerCount", C.c_uint32),
                ("ppEnabledLayerNames", C.POINTER(C.c_char_p)), ("enabledExtensionCount", C.c_uint32),
                ("ppEnabledExtensionNames", C.POINTER(C.c_char_p))]


class Extension(C.Structure):
    _fields_ = [("name", C.c_char * 256), ("version", C.c_uint32)]


class PhysicalProperties(C.Structure):
    # Vulkan 64-bit ABI prefix. The trailing uint64 storage is deliberately larger
    # than VkPhysicalDeviceLimits + VkPhysicalDeviceSparseProperties; neither is read.
    _fields_ = [("apiVersion", C.c_uint32), ("driverVersion", C.c_uint32),
                ("vendorID", C.c_uint32), ("deviceID", C.c_uint32), ("deviceType", C.c_uint32),
                ("deviceName", C.c_char * 256), ("pipelineCacheUUID", C.c_ubyte * 16),
                ("unread_tail", C.c_uint64 * 256)]


class DriverProperties(C.Structure):
    _fields_ = [("sType", C.c_uint32), ("pNext", C.c_void_p), ("driverID", C.c_uint32),
                ("driverName", C.c_char * 256), ("driverInfo", C.c_char * 256),
                ("conformance", C.c_ubyte * 4)]


class Properties2(C.Structure):
    _fields_ = [("sType", C.c_uint32), ("pNext", C.c_void_p), ("properties", PhysicalProperties)]


def loaded_images():
    # Public dyld image enumeration, not an injected library or private Metal API.
    dyld = C.CDLL(None)
    dyld._dyld_image_count.argtypes, dyld._dyld_image_count.restype = [], C.c_uint32
    dyld._dyld_get_image_name.argtypes, dyld._dyld_get_image_name.restype = [C.c_uint32], C.c_char_p
    count = dyld._dyld_image_count()
    require(0 < count < 10000, "invalid dyld image count")
    names = []
    for index in range(count):
        value = dyld._dyld_get_image_name(index)
        require(value is not None, "missing dyld image name")
        names.append(value.decode("utf-8", errors="strict"))
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loader", type=Path, required=True)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--expected-name", required=True)
    parser.add_argument("--expected-driver-info", required=True)
    args = parser.parse_args()
    require(sys.platform == "darwin" and platform.machine() == "arm64" and C.sizeof(C.c_void_p) == 8,
            "native arm64 Darwin required")
    require(PhysicalProperties.pipelineCacheUUID.offset == 276 and Properties2.properties.offset == 16,
            "Vulkan ABI prefix offsets differ")
    require(C.sizeof(PhysicalProperties) >= 2048, "property receive storage is too small")
    loader_pin, driver_pin = pin(args.loader), pin(args.driver)
    loader = C.CDLL(loader_pin["path"])
    gipa = loader.vkGetInstanceProcAddr
    gipa.argtypes, gipa.restype = [C.c_void_p, C.c_char_p], C.c_void_p

    def function(instance, name, result, parameters):
        address = gipa(instance, name.encode())
        require(address, "Vulkan entry point missing: " + name)
        return C.CFUNCTYPE(result, *parameters)(address)

    enum_ext = function(None, "vkEnumerateInstanceExtensionProperties", C.c_int32,
                        [C.c_char_p, C.POINTER(C.c_uint32), C.POINTER(Extension)])
    count = C.c_uint32()
    require(enum_ext(None, C.byref(count), None) == 0 and 0 < count.value < 1024, "extension count failed")
    extensions = (Extension * count.value)()
    expected_count = count.value
    require(enum_ext(None, C.byref(count), extensions) == 0 and count.value == expected_count,
            "extension enumeration changed or was incomplete")
    require(b"VK_KHR_portability_enumeration" in [bytes(item.name) for item in extensions],
            "portability enumeration is unavailable")
    extension_names = (C.c_char_p * 1)(b"VK_KHR_portability_enumeration")
    # Match original Sirius CreateInstance's requested Vulkan 1.3. MoltenVK uses
    # that major/minor with its compiled header patch in reported properties.
    application = ApplicationInfo(0, None, b"Sirius issue48 property query", 0, None, 0, (1 << 22) | (3 << 12))
    create_info = InstanceInfo(1, None, 1, C.pointer(application), 0, None, 1, extension_names)
    instance = C.c_void_p()
    create = function(None, "vkCreateInstance", C.c_int32,
                      [C.POINTER(InstanceInfo), C.c_void_p, C.POINTER(C.c_void_p)])
    require(create(C.byref(create_info), None, C.byref(instance)) == 0 and instance.value, "instance creation failed")
    destroy = function(instance, "vkDestroyInstance", None, [C.c_void_p, C.c_void_p])
    try:
        enumerate_devices = function(instance, "vkEnumeratePhysicalDevices", C.c_int32,
                                     [C.c_void_p, C.POINTER(C.c_uint32), C.POINTER(C.c_void_p)])
        count = C.c_uint32()
        require(enumerate_devices(instance, C.byref(count), None) == 0 and 0 < count.value < 64,
                "physical-device count failed")
        devices = (C.c_void_p * count.value)()
        expected_count = count.value
        require(enumerate_devices(instance, C.byref(count), devices) == 0 and count.value == expected_count,
                "physical-device enumeration changed or was incomplete")
        require(0 <= args.index < count.value, "selected index is absent")
        driver = DriverProperties()
        driver.sType = 1000196000
        properties = Properties2()
        properties.sType, properties.pNext = 1000059001, C.addressof(driver)
        query = function(instance, "vkGetPhysicalDeviceProperties2", None, [C.c_void_p, C.POINTER(Properties2)])
        query(devices[args.index], C.byref(properties))
        p = properties.properties
        text = lambda value: bytes(value).decode("utf-8", errors="strict")
        require(text(p.deviceName) == args.expected_name and text(driver.driverName) == "MoltenVK"
                and text(driver.driverInfo) == args.expected_driver_info and driver.driverID == 14,
                "queried physical device/driver differs")
        uuid = bytes(p.pipelineCacheUUID).hex()
        require(uuid != "00" * 16 and p.deviceType in (1, 2), "invalid cache UUID or nonphysical device")
        images = loaded_images()
        molten = {str(Path(path).resolve()) for path in images if Path(path).name.startswith("libMoltenVK")}
        vulkan = {str(Path(path).resolve()) for path in images if Path(path).name.startswith("libvulkan")}
        require(molten == {driver_pin["path"]} and vulkan == {loader_pin["path"]},
                "actual loaded Vulkan loader/driver differs or is ambiguous")
        require(pin(args.loader) == loader_pin and pin(args.driver) == driver_pin, "queried binaries changed")
        print(json.dumps({"passed": True, "device_count": count.value, "selected_index": args.index,
              "device": {"name": text(p.deviceName), "vendor_id": p.vendorID, "device_id": p.deviceID,
                         "api_version": p.apiVersion, "driver_version": p.driverVersion,
                         "device_type": p.deviceType, "driver_id": driver.driverID,
                         "driver_name": text(driver.driverName), "driver_info": text(driver.driverInfo),
                         "pipeline_cache_uuid": uuid},
              "loader": loader_pin, "driver": driver_pin, "loaded_images": images,
              "scope": "Physical-device properties and loaded images only; no VkDevice, pipeline or dispatch",
              "qualification_claimed": False}))
    finally:
        destroy(instance, None)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(json.dumps({"passed": False, "error": f"{type(error).__name__}: {error}", "qualification_claimed": False}))
        raise SystemExit(1)
