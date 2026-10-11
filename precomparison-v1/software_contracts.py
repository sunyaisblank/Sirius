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
 'timestamps': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
 'wide': 'RetainedComputeTest.WideDeviceTimestampsPreserveOriginalIntervalResults',
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

def timestamp_result(props, producer, wide=False):
    software = True
    assert producer in ('baseline','candidate')
    assert props['selected_device'] == ('llvmpipe (LLVM 20.1.2, 256 bits)' if software else 'AMD Radeon 780M Graphics')
    assert props['selected_driver'] == ('llvmpipe: Mesa 25.2.8-0ubuntu0.24.04.2 (LLVM 20.1.2)' if software else 'AMD proprietary driver: 26.8.1 (LLPC)')
    assert props['coupled_rows'] == '12' and props['coupled_projection_budget'] == '24'
    assert props['product_mode'] == ('fp64' if wide else 'default')
    assert props.get('wide_exercised') == ('1' if wide else None)
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
    # Historical shared observer metadata provides stage/complete-size/order authority;
    # current baseline/candidate digest equality is checked between all slots.
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
    arrays = document(WORK / ('baseline.json' if producer == 'baseline' else 'build-and-arrays.json'))['arrays']
    suffix = 'PortableFp64Shader' if wide else 'PortableShader'
    names = ['kCamera'+suffix, 'kTransportPortableNormalSumShader', 'kEndpointPortableNormalSumShader',
             'kDense'+suffix, 'kInitialize'+suffix, 'kRayCamera'+suffix, 'kDopriPhasePortableNormalSumShader']
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
