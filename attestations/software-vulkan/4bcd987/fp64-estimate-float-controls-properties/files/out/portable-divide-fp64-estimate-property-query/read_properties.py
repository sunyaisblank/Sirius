"""Admit one completed property observation; never loads or calls Vulkan."""
from pathlib import Path
import argparse
import hashlib
import json
import struct

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REVISION = '4498cc41168ad29fc4c111819a4cb7ccaa54a8bf'


def load(path):
    return json.loads(path.read_text())


def pin(path):
    path = Path(path)
    resolved = path.resolve(strict=True)
    digest = hashlib.sha256()
    with resolved.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'resolved_path': str(resolved), 'bytes': resolved.stat().st_size, 'sha256': digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    assert Path.cwd() == ROOT and args.output.resolve(strict=True) == HERE / 'read-properties'
    folder = HERE / 'query'
    owner = load(folder / 'owner.json')
    action = load(HERE / 'actions.json')['query']
    assert owner['source'] == {'revision': REVISION, 'status': ''}
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == []
    assert owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True
    assert owner['observed_births'] and all(row['absent'] is True for row in owner['observed_births'])
    assert owner['limits'] == {'seconds': 60, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    assert owner['argv'] == [item.replace('{folder}', str(folder)) for item in action['argv']]
    assert owner['wall_seconds'] < 60 and owner['peak_sampled_rss_kib'] <= 4194304
    assert owner['minimum_sampled_available_kib'] >= 2097152 and owner['provider_mappings']
    assert owner['declared_output_environment'] == {'GTEST_OUTPUT': '',
        'VK_DRIVER_FILES': '/usr/share/vulkan/icd.d/lvp_icd.json', 'SIRIUS_VULKAN_DEVICE': '0'}
    assert load(folder / 'inputs-before.json') == load(folder / 'inputs-after.json')
    frozen = load(HERE / 'source-freeze.json')['pins']
    for name, expected in frozen.items():
        assert pin(ROOT / name) == expected, name
    observed = load(folder / 'property-query.json')
    assert observed['schema'] == 'sirius-public-float-controls-property-query-v1'
    assert observed['source_revision'] == REVISION and observed['observation_complete'] is True
    assert observed['instance_created'] is True and observed['instance_destroyed'] is True
    assert all(observed[name] is False for name in ['logical_device_created', 'shader_pipeline_dispatch_or_WSI',
        'estimate_accuracy_established', 'logical_device_float64_enablement_observed'])
    assert [row['name'] for row in observed['api_calls']] == ['vkCreateInstance',
        'vkEnumeratePhysicalDevices', 'vkEnumeratePhysicalDevices', 'vkGetPhysicalDeviceProperties2',
        'vkGetPhysicalDeviceFeatures', 'vkDestroyInstance']
    assert all(row.get('result', 0) == 0 for row in observed['api_calls'])
    assert observed['instance_request']['api_version'] == (1 << 22) | (2 << 12)
    assert observed['instance_request']['flags'] == 0
    assert observed['instance_request']['layers'] == observed['instance_request']['extensions'] == []
    abi = load(HERE / 'abi-layout.json')
    for row in [*observed['raw_structs'].values(), observed['instance_request']['application_raw'],
                observed['instance_request']['create_info_raw']]:
        assert pin(folder / row['file']) == {key: row[key] for key in ['resolved_path', 'bytes', 'sha256']}
    raw = (folder / 'float-controls.bin').read_bytes()
    assert len(raw) == 88 and struct.unpack_from('<I', raw)[0] == 1000197000
    values = {row['name']: struct.unpack_from('<I', raw, row['offset'])[0]
        for row in abi['structs']['VkPhysicalDeviceFloatControlsProperties']['fields'][2:]}
    assert values == observed['float_controls']
    raw_features = (folder / 'core-features.bin').read_bytes()
    assert len(raw_features) == 220
    features = {row['name']: struct.unpack_from('<I', raw_features, row['offset'])[0]
        for row in abi['structs']['VkPhysicalDeviceFeatures']['fields']}
    assert features == observed['supported_core_features']
    selected = observed['selected']
    assert selected['index'] == 0 and 0 < selected['enumerated_count'] <= 64
    assert [row['count'] for row in observed['api_calls'] if row['name'] == 'vkEnumeratePhysicalDevices'] == [selected['enumerated_count']] * 2
    raw_base = (folder / 'properties2.bin').read_bytes()
    assert len(raw_base) == 840
    assert list(struct.unpack_from('<5I', raw_base, 16)) == [selected[key] for key in
        ['api_version', 'driver_version', 'vendor_id', 'device_id', 'device_type']]
    assert raw_base[36:292].split(b'\0', 1)[0].decode() == selected['name']
    assert raw_base[292:308].hex() == selected['pipeline_cache_uuid']
    raw_id = (folder / 'id-properties.bin').read_bytes()
    assert len(raw_id) == 64 and struct.unpack_from('<I', raw_id)[0] == 1000071004
    assert raw_id[16:32].hex() == selected['device_uuid'] and raw_id[32:48].hex() == selected['driver_uuid']
    assert raw_id[48:56].hex() == selected['device_luid']
    assert list(struct.unpack_from('<2I', raw_id, 56)) == [selected['device_node_mask'], selected['device_luid_valid']]
    raw_driver = (folder / 'driver-properties.bin').read_bytes()
    assert len(raw_driver) == 536 and struct.unpack_from('<I', raw_driver)[0] == 1000196000
    assert struct.unpack_from('<I', raw_driver, 16)[0] == selected['driver_id']
    assert raw_driver[20:276].split(b'\0', 1)[0].decode() == selected['driver_name']
    assert raw_driver[276:532].split(b'\0', 1)[0].decode() == selected['driver_info']
    assert list(raw_driver[532:536]) == selected['conformance_version']
    images = observed['loaded_images']
    assert images['captured_with_instance_alive'] is True
    context = images['process']
    assert context['pid'] == context['session'] == owner['child_pid']
    assert str(context['start_ticks']) == str(owner['child_start_ticks']) and context['boot_id'] == owner['boot_id']
    assert pin(folder / 'maps-while-alive.txt') == images['raw_maps']
    raw_maps = (folder / 'maps-while-alive.txt').read_text().splitlines()
    assert set(images['images']) == {'loader', 'provider', 'llvm'}
    for label, image in images['images'].items():
        assert {key: image[key] for key in ['resolved_path', 'bytes', 'sha256']} == frozen[image['resolved_path']]
        assert image['mappings'] and any('x' in row['permissions'] for row in image['mappings'])
        assert all(row['raw'] in raw_maps and row['path'] == image['resolved_path']
            and row['inode'] == image['inode']
            and tuple(int(part, 16) for part in row['dev'].split(':')) == (image['device_major'], image['device_minor'])
            for row in image['mappings'])
    checks = {'physical_api_at_least_1_2': selected['api_version'] >= (1 << 22) | (2 << 12),
        'shaderFloat64_supported': features['shaderFloat64'] == 1,
        'shaderRoundingModeRTEFloat32': values['shaderRoundingModeRTEFloat32'] == 1,
        'roundingModeIndependence_ALL_or_32_BIT_ONLY': values['roundingModeIndependence'] in (0, 1)}
    assert checks == observed['profile_checks'] and all(checks.values()) and observed['profile_admissible'] is True
    result = {'schema': 'sirius-fp64-estimate-property-admission-v1', 'passed': True,
        'source_revision': REVISION, 'profile_checks': checks, 'selected': selected,
        'float_controls': values, 'supported_core_features': features,
        'query_owner': pin(folder / 'owner.json'), 'query_result': pin(folder / 'property-query.json'),
        'inputs_before': pin(folder / 'inputs-before.json'), 'inputs_after': pin(folder / 'inputs-after.json'),
        'raw_maps': images['raw_maps'], 'source_freeze': pin(HERE / 'source-freeze.json'),
        'observed_actual_live_images': images['images'],
        'logical_device_float64_enablement_observed': False, 'estimate_accuracy_established': False,
        'future_capture_requires_exact_property_receipt_join': True,
        'qualification_claimed': False, 'scope': 'RTE32 plus internal Float64 mixed-width profile eligibility only. No RTE64 premise, logical-device enablement, arithmetic accuracy or timing result.'}
    with (args.output / 'property-accepted.json').open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps({'passed': True, 'profile_checks': checks, 'estimate_accuracy_established': False}, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
