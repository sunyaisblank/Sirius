"""Read original native controls without launching any work."""
from pathlib import Path
import hashlib, json, math, re, subprocess, xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent.parent.parent
WORK = Path(__file__).resolve().parent
CASES = {
 'dense': 'RetainedComputeTest.DenseSegmentsPreserveSmallCovariantArrivalDerivatives',
 'dense_authority': 'RetainedComputeTest.ZeroFractionArrivalsPreserveProgramAuthorityAndCompleteRefusal',
 'dopri_reference': 'RetainedDopriTest.IndependentQuarticPreservesCompletePhase',
 'dopri_refusal': 'RetainedDopriTest.MalformedRowsRefuseAndRecover',
 'dopri_projection': 'RetainedDopriTest.CurvedTransportConnectsSamplerAndProjection',
 'dopri_connected': 'RetainedDopriTest.ConnectedIntervalsPreserveTrialsAndIndependentBudgets',
 'dopri_arrival': 'RetainedDopriTest.SamplerPreservesPhysicalArrivalAndRejectsInconsistentRates',
 'timestamps': 'RetainedComputeTest.WideDeviceTimestampsPreserveOriginalIntervalResults',
}

MODERN = [('transport', 12), ('endpoint', 24), ('transport', 12), ('endpoint', 24),
          ('transport', 12), ('endpoint', 24), ('dopri_phase', 24), ('dopri_phase', 24),
          ('endpoint', 12), ('dopri_phase', 12), ('endpoint', 12)]
LEGACY = ['transport', 'endpoint', 'transport', 'endpoint', 'transport', 'endpoint+dense']
MARKER_MS = ['device_ms', 'host_submit_ms', 'host_wait_ms', 'host_completion_ms',
             'pipeline_ms', 'setup_ms', 'cleanup_ms', 'total_ms']

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def identity(path):
    return {'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': sha(path)}

def read(path):
    assert path.is_file() and path.stat().st_size <= 16 * 1024 * 1024, str(path)
    return path.read_text(encoding='utf-8-sig')

def document(path):
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    return json.loads(read(path), parse_constant=invalid)

def number(value):
    assert not isinstance(value, bool)
    result = float(value)
    assert math.isfinite(result) and result >= 0
    return result

def fields(value):
    result = {}
    for item in value.split(';'):
        key, sep, word = item.partition('=')
        assert sep and key not in result
        result[key] = word
    return result

def xml_properties(path, case):
    xml = ET.fromstring(read(path))
    assert int(xml.get('tests', '-1')) == 1
    for node in [xml, *xml.findall('.//testsuite')]:
        for key in ('failures', 'errors', 'skipped', 'disabled'):
            assert int(node.get(key, '0')) == 0
    assert not any(xml.findall('.//' + tag) for tag in ('failure', 'error', 'skipped'))
    cases = xml.findall('.//testcase')
    assert len(cases) == 1 and cases[0].get('classname') + '.' + cases[0].get('name') == CASES[case]
    assert cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
    props = {}
    for item in xml.findall('.//property'):
        key = item.get('name')
        assert key and key not in props
        props[key] = item.get('value')
    return props

def timestamp_result(props, producer):
    software = producer == 'software'
    assert props['selected_device'] == ('llvmpipe (LLVM 20.1.2, 256 bits)' if software else 'AMD Radeon 780M Graphics')
    assert props['selected_driver'] == ('llvmpipe: Mesa 25.2.8-0ubuntu0.24.04.2 (LLVM 20.1.2)' if software else 'AMD proprietary driver: 26.8.1 (LLPC)')
    assert props['coupled_rows'] == '12' and props['coupled_projection_budget'] == '24'
    assert props['product_mode'] == 'fp64' and props['wide_exercised'] == '1'
    assert props['retained_scalar_route'] == ('portable with guarded normal sums' if software else 'native binary32')
    for key, expected in [('fma_fp32_enabled', not software), ('preserves_fp32_signed_zero_inf_nan', True)]:
        assert props[key] in ('0', '1', 'false', 'true')
        assert (props[key] in ('1', 'true')) is expected
    assert props['device_observations'] == '51'
    bits = int(props['timestamp_valid_bits']); period = number(props['timestamp_period_ns'])
    assert bits == 64 and period == (1 if software else 10)
    expected = {f'device_observation_{i}' for i in range(51)}
    assert {key for key in props if key.startswith('device_observation_')} == expected
    order = [(f'legacy_repeat{r}', stage) for r in range(3) for stage in LEGACY]
    order += [(f'coupled_repeat{r}', f'{stage}:rows={rows}') for r in range(3) for stage, rows in MODERN]
    phases = {}
    for index, pair in enumerate(order):
        marker = fields(props[f'device_observation_{index}'])
        assert set(marker) == {'phase', 'stage', 'begin', 'end', 'available0', 'available1', *MARKER_MS}
        assert (marker['phase'], marker['stage']) == pair
        assert marker['available0'] == marker['available1'] == '1'
        begin, end = int(marker['begin']), int(marker['end'])
        assert 0 <= begin < 2**64 and 0 <= end < 2**64
        for key in MARKER_MS:
            marker[key] = number(marker[key])
        expected_ms = ((end - begin) & ((1 << bits) - 1)) * period / 1e6
        assert math.isclose(marker['device_ms'], expected_ms, rel_tol=1e-9, abs_tol=1e-9)
        assert math.isclose(marker['host_completion_ms'], marker['host_submit_ms'] + marker['host_wait_ms'], rel_tol=1e-9, abs_tol=1e-9)
        assert math.isclose(marker['total_ms'], marker['pipeline_ms'] + marker['setup_ms'] + marker['host_completion_ms'] + marker['cleanup_ms'], rel_tol=1e-9, abs_tol=1e-9)
        phases.setdefault(pair[0], []).append(marker)
    readbacks = {}
    expected_reads = {f'legacy_repeat{r}_readback_{i}' for r in range(3) for i in range(7)}
    expected_reads |= {f'coupled_readback_{i}' for i in range(11)}
    actual_reads = {k for k in props if re.fullmatch(r'legacy_repeat\d+_readback_\d+|coupled_readback_\d+', k)}
    assert actual_reads == expected_reads
    # Original software XML provides only stage/complete-size/order authority;
    # native digest equality is checked between native producers separately.
    reference = document(WORK / 'expected-readback-layout.json')['metadata']
    for key in expected_reads:
        value = fields(props[key]); original = reference[key]
        assert set(value) == {'stage', 'bytes', 'sha256'}
        assert value['stage'] == original['stage'] and value['bytes'] == original['bytes']
        assert re.fullmatch('[0-9a-f]{64}', value['sha256'])
        readbacks[key] = value
    for repeat in (1, 2):
        for i in range(7):
            assert readbacks[f'legacy_repeat{repeat}_readback_{i}'] == readbacks[f'legacy_repeat0_readback_{i}']
    legacy_calls = [number(props[f'legacy_repeat{r}_complete_call_ms']) for r in range(3)]
    modern_calls = [number(props[f'coupled_repeat{r}_complete_call_ms']) for r in range(3)]
    for r in range(3):
        assert math.isclose(modern_calls[r], number(props[f'coupled_repeat{r}_attempt_ms']) + number(props[f'coupled_repeat{r}_sample_ms']), rel_tol=1e-12, abs_tol=1e-12)
    result = {}
    for family, prefix, calls in [('legacy', 'legacy_repeat', legacy_calls), ('modern', 'coupled_repeat', modern_calls)]:
        records = []
        for r in range(3):
            markers = phases[prefix + str(r)]
            records.append({'complete_call_ms': calls[r],
                            'device_ms': math.fsum(m['device_ms'] for m in markers),
                            'host_completion_ms': math.fsum(m['host_completion_ms'] for m in markers),
                            'markers': markers})
        result[family] = {'repeats': records, 'sums': {k: math.fsum(v[k] for v in records) for k in ('device_ms', 'host_completion_ms', 'complete_call_ms')}}
    result['readbacks'] = readbacks
    result['identity'] = {key: props[key] for key in ('selected_device', 'selected_driver', 'retained_scalar_route', 'product_mode', 'timestamp_valid_bits', 'timestamp_period_ns', 'coupled_rows', 'coupled_projection_budget', 'timestamp_scope', 'legacy_call_timing_scope', 'coupled_call_timing_scope')}
    stages = ['camera', 'transport', 'endpoint', 'dense', 'initialize', 'ray_camera', 'dopri_phase']
    arrays = document(WORK / 'build-and-arrays.json')['arrays']
    names = ['kCameraFp64Shader', 'kTransportFmaShader', 'kEndpointFmaShader',
             'kDenseFp64Shader' if producer == 'baseline' else 'kDenseFmaShader',
             'kInitializeFp64Shader', 'kRayCameraFp64Shader',
             'kDopriPhaseFp64Shader' if producer == 'baseline' else 'kDopriPhaseFmaShader']
    if software:
        names = ['kCameraPortableFp64Shader', 'kTransportPortableNormalSumShader', 'kEndpointPortableNormalSumShader',
                 'kDensePortableFp64Shader', 'kInitializePortableFp64Shader', 'kRayCameraPortableFp64Shader',
                 'kDopriPhasePortableNormalSumShader']
    assert {k for k in props if k.startswith('selected_module_')} == {f'selected_module_{i}' for i in range(7)}
    selected = {}
    for i, (stage, name) in enumerate(zip(stages, names)):
        module = fields(props[f'selected_module_{i}'])
        assert set(module) == {'stage', 'bytes', 'sha256'}
        assert module == {'stage': stage, 'bytes': str(arrays[name]['bytes']), 'sha256': arrays[name]['sha256']}, (stage, module)
        selected[stage] = {'array': name, **module}
    result['selected_modules'] = selected
    result['expected_complete_readback_count'] = 32
    return result

def read_case(mode, label, case, bindings):
    assert mode in ('controls', 'matched') and case in CASES
    assert (mode == 'controls' and label == 'candidate') or (mode == 'matched' and case == 'timestamps' and label in ('a1', 'b1', 'b2', 'a2'))
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() == bindings['candidate_revision']
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT)
    for relative, seal in bindings['files'].items():
        path = Path(relative)
        if not path.is_absolute(): path = ROOT / path
        assert path.stat().st_size == seal['bytes'] and sha(path) == seal['sha256'], relative
    producer = 'baseline' if mode == 'baseline' or (mode == 'matched' and label.startswith('a')) else 'candidate'
    folder = WORK / 'native' / mode / label / case
    owner = document(folder / 'owner.json')
    bridge_path = WORK / ('native-' + mode + '-' + label + '-' + case + '-bridge.json')
    bridge = document(bridge_path)
    assert bridge['completed'] and bridge['accepted_owner_execution'] and bridge['returncode'] == 0
    binding_path = WORK / ('baseline-native-bindings.json' if mode == 'baseline' else 'execution-bindings.json')
    assert bridge['execution_bindings_sha256'] == sha(binding_path)
    assert bridge['input_bindings'] == bindings['files']
    assert bridge['producer'] == producer and bridge['case'] == case and bridge['mode'] == mode and bridge['label'] == label
    assert bridge['input_seals_current_and_exact_at_end'] and bridge['execution_bindings_unchanged_at_end']
    assert bridge['source_tree_clean_at_start'] and bridge['source_tree_clean_at_end']
    assert bridge['live_source_revision'] == bridge['live_source_revision_at_end'] == bindings['candidate_revision']
    assert owner['status'] == 'completed' and owner['returncode'] == 0 and not owner['outer_stop']
    assert not owner['cleanup_errors']
    for key in ('owned_process_absent', 'numerical_control_passed', 'test_executable_unchanged', 'loaded_driver_modules_recorded', 'source_build_gate_unchanged'):
        assert owner[key] is True
    assert owner['expected_cases'] == [CASES[case]]
    assert owner['outer_guard_seconds'] == 300 and owner['rss_guard_bytes'] == 4294967296
    assert owner['source_revision'] == bridge['source_revision'] == bindings['producers'][producer]['source_revision']
    assert owner['live_source_revision'] == bindings['candidate_revision']
    assert owner['source_build_gate_sha256'] == bindings['producers'][producer]['gate_sha256']
    exe = bindings['producers'][producer]['stage'] + '/tests/backend/Release/sirius_backend_tests.exe'
    assert owner['test_executable_sha256'] == bindings['files'][exe]['sha256']
    for key, word in owner['environment'].items():
        if key not in ('TEMP', 'TMP'):
            assert word in (None, '')
    for key in ('owner_pid', 'test_pid'):
        assert type(owner[key]) is int and owner[key] > 0
    for key in ('owner_start_time_utc', 'test_start_time_utc'):
        assert owner[key]
    modules = document(folder / 'driver-modules.json')
    assert isinstance(modules, list) and len(modules) == 2
    module_identities = []
    bound_paths = {path.casefold(): path for path in bindings['files'] if path.startswith('/mnt/c/')}
    assert {module['name'].lower() for module in modules} == {'vulkan-1.dll', 'amdvlk64.dll'}
    for module in modules:
        path = '/mnt/' + module['path'][0].lower() + module['path'][2:].replace('\\', '/')
        assert {'bytes': module['bytes'], 'sha256': module['sha256']} == bindings['files'][bound_paths[path.casefold()]]
        module_identities.append({k: module[k].casefold() if k in ('name', 'path') else module[k] for k in ('name', 'path', 'bytes', 'sha256', 'version')})
    props = xml_properties(folder / 'gtest.xml', case)
    record = {'case': case, 'label': label, 'source_revision': owner['source_revision'],
              'owner': owner, 'modules': sorted(module_identities, key=lambda x: x['name'].lower()),
              'artifacts': [identity(folder / name) for name in ('owner.json', 'gtest.xml', 'stdout.log', 'stderr.log', 'driver-modules.json', 'memory.json')] + [identity(bridge_path)]}
    if case == 'timestamps':
        record['timestamps'] = timestamp_result(props, producer)
    if case == 'dense':
        assert props['dense_wide_exercised'] == '1'
        assert props['dense_wide_shader_sha256'] == document(WORK / 'build-and-arrays.json')['arrays']['kDenseFmaShader']['sha256']
    return record
