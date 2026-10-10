"""Validate one original FP64 science case and its finite observation contract."""
from pathlib import Path
import hashlib
import json
import math
import re
import sys
import xml.etree.ElementTree as ET

work = Path(__file__).resolve().parent
folder = work / sys.argv[1]
native = len(sys.argv) == 3 and sys.argv[2] == 'native'
assert len(sys.argv) == (3 if native else 2)
source = json.loads((work / 'source.json').read_text())['source_revision']
expected = json.loads((work / 'expected-stream.json').read_text())
arrays = json.loads((work / 'build-and-arrays.json').read_text())['arrays']
owner = json.loads((folder / 'owner.json').read_text(encoding='utf-8-sig'))
if native:
    assert owner['status'] == 'completed' and owner['returncode'] == 0 and not owner['outer_stop']
    assert owner['source_revision'] == owner['live_source_revision'] == source
    assert owner['owned_process_absent'] and owner['numerical_control_passed']
    assert owner['test_executable_unchanged'] and owner['source_build_gate_unchanged']
    assert owner['loaded_driver_modules_recorded'] and not owner['cleanup_errors']
else:
    assert owner['passed'] and owner['child_exit'] == 0 and owner['stop_reason'] is None
    assert owner['source'] == {'revision': source, 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['remaining_owned_processes'] == [] and not owner['cleanup_errors']
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
for file in (('gtest.xml',) if native else ('gtest.xml', 'ctest.xml')):
    tree = ET.parse(folder / file)
    cases = tree.findall('.//testcase')
    assert len(cases) == 1
    assert not tree.findall('.//failure') and not tree.findall('.//error') and not tree.findall('.//skipped')
    for element in tree.iter():
        for key in ('failures', 'errors', 'skipped', 'disabled'):
            if key in element.attrib: assert int(element.attrib[key]) == 0
    case = cases[0]
    assert case.attrib['name'] in ('Fp64ProductsPreserveIndependentScienceOrDeclineUnsupportedDevices',
        'RetainedComputeTest.Fp64ProductsPreserveIndependentScienceOrDeclineUnsupportedDevices')
    if file == 'gtest.xml':
        assert case.attrib['status'] == 'run' and case.attrib['result'] == 'completed'
        properties = case.findall('./properties/property')
        assert len({p.attrib['name'] for p in properties}) == len(properties)
        props = {p.attrib['name']: p.attrib['value'] for p in properties}
        gtest_seconds = float(case.attrib['time'])

def fields(value):
    pairs = [item.split('=', 1) for item in value.split(';')]
    assert len({key for key, _ in pairs}) == len(pairs)
    return dict(pairs)

def flag(key):
    value = props[key]
    assert value in ('0', '1', 'false', 'true'), (key, value)
    return value in ('1', 'true')

if native:
    plan = json.loads((work / 'native-plan.json').read_text())
    assert plan['source_revision'] == source
    for key in ('selected_device', 'selected_driver'):
        assert props[key] == plan[key], (key, props[key])
    for key, value in plan['required_capabilities'].items():
        assert flag(key) is value, (key, props[key])

assert flag('supports_fp64') and flag('rounds_fp64_to_nearest'), 'refusal is not supported-path observation'
assert props['timestamp_observation'] == 'enabled', 'science may pass without timing; observation acceptance may not'
assert int(props['timestamp_valid_bits']) in range(36, 65) and float(props['timestamp_period_ns']) > 0
modules = []
portable = not flag('preserves_fp32_denormals') or not flag('rounds_fp32_to_nearest')
normal_sum = portable and flag('rounds_fp32_to_nearest')
fma = not portable and flag('fma_fp32_enabled')
base_names = ['Camera', 'Transport', 'Endpoint', 'Dense', 'Initialize', 'RayCamera', 'DopriPhase']
for index, (stage, base) in enumerate(zip(expected['module_order'], base_names)):
    actual = fields(props[f'selected_module_{index}'])
    assert actual['stage'] == stage
    name = 'k' + base + ('PortableNormalSum' if normal_sum and base in ('Transport', 'Endpoint', 'DopriPhase')
                        else 'Fma' if fma and base in ('Transport', 'Endpoint')
                        else 'PortableFp64' if portable else 'Fp64') + 'Shader'
    assert {'bytes': int(actual['bytes']), 'sha256': actual['sha256']} == arrays[name], (stage, name, actual)
    modules.append({'stage': stage, 'array': name, **arrays[name]})
readbacks = []
call_ms = {}
for phase, stream in expected['readbacks'].items():
    assert int(props[phase + '_readbacks']) == len(stream)
    observed = []
    for index, pair in enumerate(stream):
        actual = fields(props[f'{phase}_readback_{index}'])
        assert [actual['stage'], int(actual['rows'])] == pair
        words = {'ray_camera': 3216, 'transport': 4699, 'endpoint': 3829, 'dense': 4279, 'initialize': 2014}
        assert int(actual['bytes']) == int(actual['rows']) * 4 * words[actual['stage']]
        assert re.fullmatch('[0-9a-f]{64}', actual['sha256'])
        observed.append(actual)
        readbacks.append({'phase': phase, 'index': index, **actual})
    call_ms[phase] = float(props[phase + '_complete_call_ms'])
    assert math.isfinite(call_ms[phase]) and call_ms[phase] >= 0
assert len(readbacks) == expected['full_active_readback_count']
assert int(props['device_observations']) == expected['wide_marker_count']
marker_stream = []
for phase, stream in expected['readbacks'].items():
    if phase == 'paired':
        marker_stream += [(phase, stage, rows) for stage, rows in stream[:-2]]
        marker_stream.append((phase, 'endpoint+dense', (4, 2)))
    else:
        marker_stream += [(phase, stage, rows) for stage, rows in stream]
markers = []
for index, (phase, stage, rows) in enumerate(marker_stream):
    actual = fields(props[f'device_observation_{index}'])
    assert (actual['phase'], actual['stage']) == (phase, stage)
    if stage == 'endpoint+dense':
        assert (int(actual['endpoint_rows']), int(actual['dense_rows'])) == rows
    else:
        # The frozen generated ABI has one row per workgroup.
        assert int(actual['groups_x']) == rows
    assert int(actual['query_result']) == 0
    assert int(actual['available0']) != 0 and int(actual['available1']) != 0
    for key in ('pipeline_ms', 'setup_ms', 'host_completion_ms', 'cleanup_ms', 'total_ms',
                'host_submit_ms', 'host_wait_ms', 'device_ms'):
        value = float(actual[key]); assert math.isfinite(value) and value >= 0
    assert math.isclose(float(actual['host_submit_ms']) + float(actual['host_wait_ms']),
                        float(actual['host_completion_ms']), rel_tol=1e-9, abs_tol=1e-9)
    mask = (1 << int(props['timestamp_valid_bits'])) - 1
    span = ((int(actual['end']) - int(actual['begin'])) & mask) * float(props['timestamp_period_ns']) * 1e-6
    assert math.isclose(span, float(actual['device_ms']), rel_tol=1e-12, abs_tol=1e-9)
    markers.append(actual)
assert len(markers) == 30
summary = {
    'source_revision': source, 'case': 'RetainedComputeTest.Fp64ProductsPreserveIndependentScienceOrDeclineUnsupportedDevices',
    'passed': True, 'native_original_body_execution': native, 'ctest_wrapper_execution': not native, 'gtest_seconds': gtest_seconds, 'owner_seconds': owner['elapsed_seconds'] if native else owner['wall_seconds'],
    'device': props['selected_device'], 'driver': props['selected_driver'],
    'selected_modules': modules, 'ordered_wide_markers': markers,
    'complete_active_readbacks': readbacks, 'call_ms': call_ms,
    'observation_scope': props['observation_scope'],
    'scope': 'One original supported-path finite scientific pass and complete wide-only observations; narrow host span only. No frame, modern-DP, coldness, isolated arithmetic/driver, speedup, full-estate or release qualification. Native acceptance also requires separate exact producer/provider/binding and terminal audits.'}
(folder / 'accepted-observations.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({k: summary[k] for k in ('source_revision', 'passed', 'gtest_seconds', 'owner_seconds', 'device', 'driver', 'call_ms')}, indent=2))
