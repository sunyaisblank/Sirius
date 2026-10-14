"""One public Vulkan1.2 physical-property query; no logical device or work."""
from pathlib import Path
import argparse
import ctypes as C
import hashlib
import json
import os
import platform
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LOADER = Path('/usr/lib/x86_64-linux-gnu/libvulkan.so.1.3.275')
PROVIDER = Path('/usr/lib/x86_64-linux-gnu/libvulkan_lvp.so')
LLVM = Path('/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1')
API12 = (1 << 22) | (2 << 12)


def require(value, message):
    if not value:
        raise RuntimeError(message)


def pin(path):
    path = Path(path)
    resolved = path.resolve(strict=True)
    before = resolved.stat()
    digest = hashlib.sha256()
    with resolved.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    after = resolved.stat()
    require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) ==
            (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), 'input changed while hashing')
    return {'resolved_path': str(resolved), 'bytes': after.st_size, 'sha256': digest.hexdigest()}


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def process_context():
    raw = Path('/proc/self/stat').read_text()
    fields = raw[raw.rindex(')') + 2:].split()
    return {'pid': os.getpid(), 'ppid': int(fields[1]), 'pgrp': int(fields[2]),
            'session': int(fields[3]), 'start_ticks': fields[19],
            'boot_id': Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
            'stat_raw': raw}


def layouts():
    require(sys.platform == 'linux' and platform.machine() == 'x86_64' and sys.byteorder == 'little',
            'exact Linux x86_64 little-endian ABI required')
    require(C.sizeof(C.c_void_p) == C.sizeof(C.c_size_t) == 8 and C.sizeof(C.c_int) == 4,
            'LP64/default C enum ABI differs')
    abi = json.loads((HERE / 'abi-layout.json').read_text())
    names = ['VkApplicationInfo', 'VkInstanceCreateInfo', 'VkPhysicalDeviceFeatures',
        'VkPhysicalDeviceLimits', 'VkPhysicalDeviceSparseProperties', 'VkPhysicalDeviceProperties',
        'VkPhysicalDeviceProperties2', 'VkPhysicalDeviceIDProperties',
        'VkConformanceVersion', 'VkPhysicalDeviceDriverProperties', 'VkPhysicalDeviceFloatControlsProperties']
    require(set(abi['structs']) == set(names), 'fixed struct schema differs')
    types = {'uint32_t': C.c_uint32, 'int32_t': C.c_int32, 'float': C.c_float,
        'VkDeviceSize': C.c_uint64, 'size_t': C.c_size_t, 'VkSampleCountFlags': C.c_uint32,
        'VkBool32': C.c_uint32, 'VkStructureType': C.c_uint32, 'VkPhysicalDeviceType': C.c_uint32,
        'VkDriverId': C.c_uint32, 'VkShaderFloatControlsIndependence': C.c_uint32,
        'VkInstanceCreateFlags': C.c_uint32, 'char': C.c_char, 'uint8_t': C.c_uint8,
        'void*': C.c_void_p, 'const void*': C.c_void_p, 'const char*': C.c_char_p,
        'const char* const*': C.POINTER(C.c_char_p)}
    for name in names:
        schema = abi['structs'][name]
        fields = [(row['name'], types[row['c_type']] * row['elements'] if row['elements'] > 1
                   else types[row['c_type']]) for row in schema['fields']]
        kind = type(name, (C.Structure,), {'_fields_': fields})
        require(C.sizeof(kind) == schema['bytes'] and C.alignment(kind) == schema['alignment'],
                'ctypes struct size/alignment differs: ' + name)
        require(all(getattr(kind, row['name']).offset == row['offset'] for row in schema['fields']),
                'ctypes field offsets differ: ' + name)
        types[name] = kind
        if name == 'VkApplicationInfo':
            types['const VkApplicationInfo*'] = C.POINTER(kind)
    return abi, types


def raw_struct(folder, name, value):
    raw = C.string_at(C.addressof(value), C.sizeof(value))
    path = folder / (name + '.bin')
    with path.open('xb') as stream:
        stream.write(raw)
    return {'file': path.name, 'address': C.addressof(value), **pin(path)}


def mapped_images(folder, expected):
    context_before = process_context()
    raw = Path('/proc/self/maps').read_text()
    (folder / 'maps-while-alive.txt').write_text(raw)
    parsed = []
    for line in raw.splitlines():
        fields = line.split(None, 5)
        if len(fields) == 6 and fields[5].startswith('/'):
            parsed.append({'range': fields[0], 'permissions': fields[1], 'offset': fields[2],
                'dev': fields[3], 'inode': int(fields[4]), 'path': fields[5], 'raw': line})
    result = {}
    for label, path in [('loader', LOADER), ('provider', PROVIDER), ('llvm', LLVM)]:
        rows = [row for row in parsed if row['path'] == str(path)]
        require(rows and any('x' in row['permissions'] for row in rows), 'actual mapped image absent: ' + label)
        info = path.stat()
        device = (os.major(info.st_dev), os.minor(info.st_dev))
        require(all(tuple(int(part, 16) for part in row['dev'].split(':')) == device
                    and row['inode'] == info.st_ino for row in rows), 'mapped path/dev/inode differs from file: ' + label)
        actual = pin(path)
        require(actual == expected[str(path)], 'actual mapped image differs from frozen provider: ' + label)
        result[label] = {**actual, 'device_major': device[0], 'device_minor': device[1],
            'inode': info.st_ino, 'mappings': rows}
    current = Path('/proc/self/maps').read_text().splitlines()
    require(all(row['raw'] in current for image in result.values() for row in image['mappings']),
            'an observed image mapping changed while hashing')
    context_after = process_context()
    require(all(context_before[key] == context_after[key] for key in ['pid', 'start_ticks', 'session', 'boot_id']),
            'query process identity changed')
    return {'process': context_before, 'images': result,
        'raw_maps': pin(folder / 'maps-while-alive.txt'), 'captured_with_instance_alive': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    folder = args.output.resolve(strict=True)
    require(Path.cwd() == ROOT and folder == HERE / 'query', 'owned query output folder required')
    require(os.environ.get('VK_DRIVER_FILES') == '/usr/share/vulkan/icd.d/lvp_icd.json'
            and os.environ.get('SIRIUS_VULKAN_DEVICE') == '0'
            and not os.environ.get('MESA_SHADER_CACHE_DISABLE'), 'original default cache/LVP/device0 environment required')
    frozen = json.loads((HERE / 'source-freeze.json').read_text())['pins']
    for name, expected in frozen.items():
        require(pin(ROOT / name) == expected, 'frozen property-query input changed: ' + name)
    abi, types = layouts()
    report = {'schema': 'sirius-public-float-controls-property-query-v1',
        'source_revision': '4498cc41168ad29fc4c111819a4cb7ccaa54a8bf',
        'observation_complete': False, 'profile_admissible': False, 'instance_created': False,
        'instance_destroyed': False, 'logical_device_created': False,
        'shader_pipeline_dispatch_or_WSI': False, 'estimate_accuracy_established': False,
        'logical_device_float64_enablement_observed': False, 'api_calls': [],
        'environment': {name: os.environ.get(name) for name in ['VK_DRIVER_FILES', 'SIRIUS_VULKAN_DEVICE', 'MESA_SHADER_CACHE_DISABLE']}}
    error = None
    try:
        loader = C.CDLL(str(LOADER))
        create, destroy = loader.vkCreateInstance, loader.vkDestroyInstance
        create.argtypes = [C.POINTER(types['VkInstanceCreateInfo']), C.c_void_p, C.POINTER(C.c_void_p)]
        create.restype = C.c_int32
        destroy.argtypes, destroy.restype = [C.c_void_p, C.c_void_p], None
        enumerate_devices = loader.vkEnumeratePhysicalDevices
        enumerate_devices.argtypes = [C.c_void_p, C.POINTER(C.c_uint32), C.POINTER(C.c_void_p)]
        enumerate_devices.restype = C.c_int32
        properties2 = loader.vkGetPhysicalDeviceProperties2
        properties2.argtypes, properties2.restype = [C.c_void_p, C.POINTER(types['VkPhysicalDeviceProperties2'])], None
        features_query = loader.vkGetPhysicalDeviceFeatures
        features_query.argtypes, features_query.restype = [C.c_void_p, C.POINTER(types['VkPhysicalDeviceFeatures'])], None
        application = types['VkApplicationInfo'](0, None, b'Sirius bounded FloatControls property query', 0, None, 0, API12)
        create_info = types['VkInstanceCreateInfo'](1, None, 0, C.pointer(application), 0, None, 0, None)
        report['instance_request'] = {'api_version': API12, 'flags': 0, 'layers': [], 'extensions': [],
            'application_raw': raw_struct(folder, 'application-info', application),
            'create_info_raw': raw_struct(folder, 'instance-create-info', create_info)}
        instance = C.c_void_p()
        result = create(C.byref(create_info), None, C.byref(instance))
        report['api_calls'].append({'name': 'vkCreateInstance', 'result': result})
        require(result == 0 and instance.value, 'public Vulkan1.2 instance creation failed')
        report['instance_created'] = True
        try:
            count = C.c_uint32()
            result = enumerate_devices(instance, C.byref(count), None)
            report['api_calls'].append({'name': 'vkEnumeratePhysicalDevices', 'result': result, 'count': count.value})
            require(result == 0 and 0 < count.value <= 64, 'physical-device count failed or exceeds bound')
            expected_count = count.value
            devices = (C.c_void_p * expected_count)()
            result = enumerate_devices(instance, C.byref(count), devices)
            report['api_calls'].append({'name': 'vkEnumeratePhysicalDevices', 'result': result, 'count': count.value})
            require(result == 0 and count.value == expected_count, 'enumeration incomplete or changed')
            physical = devices[0]
            require(physical, 'selected device0 handle absent')
            driver = types['VkPhysicalDeviceDriverProperties']()
            driver.sType = abi['sTypes']['VkPhysicalDeviceDriverProperties']
            identifiers = types['VkPhysicalDeviceIDProperties']()
            identifiers.sType, identifiers.pNext = abi['sTypes']['VkPhysicalDeviceIDProperties'], C.addressof(driver)
            controls = types['VkPhysicalDeviceFloatControlsProperties']()
            controls.sType, controls.pNext = abi['sTypes']['VkPhysicalDeviceFloatControlsProperties'], C.addressof(identifiers)
            properties = types['VkPhysicalDeviceProperties2']()
            properties.sType, properties.pNext = abi['sTypes']['VkPhysicalDeviceProperties2'], C.addressof(controls)
            properties2(physical, C.byref(properties))
            report['api_calls'].append({'name': 'vkGetPhysicalDeviceProperties2'})
            features = types['VkPhysicalDeviceFeatures']()
            features_query(physical, C.byref(features))
            report['api_calls'].append({'name': 'vkGetPhysicalDeviceFeatures'})
            report['raw_structs'] = {name: raw_struct(folder, name, value) for name, value in
                [('properties2', properties), ('float-controls', controls), ('id-properties', identifiers),
                 ('driver-properties', driver), ('core-features', features)]}
            require(properties.sType == 1000059001 and properties.pNext == C.addressof(controls)
                    and controls.sType == 1000197000 and controls.pNext == C.addressof(identifiers)
                    and identifiers.sType == 1000071004 and identifiers.pNext == C.addressof(driver)
                    and driver.sType == 1000196000 and not driver.pNext, 'property chain was changed')
            p = properties.properties
            text = lambda value: bytes(value).decode('utf-8', errors='strict')
            report['selected'] = {'index': 0, 'enumerated_count': count.value, 'physical_handle': physical,
                'name': text(p.deviceName), 'vendor_id': p.vendorID, 'device_id': p.deviceID,
                'device_type': p.deviceType, 'api_version': p.apiVersion, 'driver_version': p.driverVersion,
                'pipeline_cache_uuid': bytes(p.pipelineCacheUUID).hex(),
                'device_uuid': bytes(identifiers.deviceUUID).hex(), 'driver_uuid': bytes(identifiers.driverUUID).hex(),
                'device_luid': bytes(identifiers.deviceLUID).hex(), 'device_luid_valid': identifiers.deviceLUIDValid,
                'device_node_mask': identifiers.deviceNodeMask, 'driver_id': driver.driverID,
                'driver_name': text(driver.driverName), 'driver_info': text(driver.driverInfo),
                'conformance_version': [getattr(driver.conformanceVersion, name)
                    for name in ['major', 'minor', 'subminor', 'patch']]}
            feature_names = [row['name'] for row in abi['structs']['VkPhysicalDeviceFeatures']['fields']]
            report['supported_core_features'] = {name: getattr(features, name) for name in feature_names}
            control_names = [row['name'] for row in abi['structs']['VkPhysicalDeviceFloatControlsProperties']['fields']][2:]
            report['float_controls'] = {name: getattr(controls, name) for name in control_names}
            require(all(value in (0, 1) for value in report['supported_core_features'].values()), 'invalid VkBool32 feature')
            require(all(value in (0, 1) for name, value in report['float_controls'].items() if name.startswith('shader')),
                    'invalid VkBool32 FloatControls')
            require(controls.roundingModeIndependence in (0, 1, 2) and controls.denormBehaviorIndependence in (0, 1, 2),
                    'unknown FloatControls independence enum')
            report['loaded_images'] = mapped_images(folder, frozen)
            report['profile_checks'] = {'physical_api_at_least_1_2': p.apiVersion >= API12,
                'shaderFloat64_supported': features.shaderFloat64 == 1,
                'shaderRoundingModeRTEFloat32': controls.shaderRoundingModeRTEFloat32 == 1,
                'roundingModeIndependence_ALL_or_32_BIT_ONLY': controls.roundingModeIndependence in (0, 1)}
            report['profile_admissible'] = all(report['profile_checks'].values())
            report['observation_complete'] = True
        finally:
            destroy(instance, None)
            report['api_calls'].append({'name': 'vkDestroyInstance'})
            report['instance_destroyed'] = True
    except BaseException as caught:
        error = caught
        report['error'] = repr(caught)
        report['observation_complete'] = False
        report['profile_admissible'] = False
    report['claim_limits'] = ['Physical-device property/support query only; no logical-device enablement, shader validation or numerical execution.',
        'RTE64 and denorm properties are preserved observations, not admission requirements for this candidate.',
        'Mixed-width profile validity does not establish division-estimate accuracy or original arithmetic acceptance.',
        'Future capture must separately bind these exact result/owner/inputs/maps and actual shaderFloat64-enabled device arguments.']
    write(folder / 'property-query.json', report)
    print(json.dumps({'observation_complete': report['observation_complete'],
        'profile_admissible': report['profile_admissible'], 'estimate_accuracy_established': False}, sort_keys=True), flush=True)
    if error is not None:
        raise error


if __name__ == '__main__':
    main()
