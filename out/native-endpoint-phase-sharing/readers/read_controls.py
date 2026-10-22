"""Admit the already completed original controls; never execute a workload."""
from pathlib import Path
from types import SimpleNamespace
import hashlib, importlib.util, json, re
R = Path.cwd(); W = Path(__file__).resolve().parent; RUN = W.parent / 'runtime'
spec = importlib.util.spec_from_file_location('original_controls', W / 'native_evidence.py')
old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
actions = doc(RUN / 'actions.json')
provider = doc(W.parent / 'providers.json')['providers']['candidate']
control8 = doc(W / 'control-8-acceptance.json'); benefit = doc(W / 'benefit.json')
decimal = doc(W / 'decimal-agreement.json')
assert control8['passed'] and benefit['all22_numeric_pass']
assert decimal['agrees_on_all22'] and decimal['all22_pass']
modules = doc(W.parent / 'compiler-result.json')['modules']
native = dict(zip(old.STAGES, ['k' + p + 'Shader' for p in old.PREFIXES]))
native.update(endpoint='kEndpointNativeOwnedShader', dopri_phase='kDopriPhaseNativeOwnedShader')
arrays = {}
C = R / 'bin/linux-gcc/src/sirius/backend/retained'
for stage, prefix in zip(old.STAGES, old.PREFIXES):
    for suffix, stem in [('Shader', ''), ('Fp64Shader', '_fp64')]:
        arrays['k' + prefix + suffix] = pin(C / ('retained_' + stage + stem + '.spv'))
for stage, prefix in [('transport', 'Transport'), ('endpoint', 'Endpoint')]:
    arrays['k' + prefix + 'FmaShader'] = pin(C / ('retained_' + stage + '_fma.spv'))
for kind in ['Endpoint', 'DopriPhase']:
    arrays['k' + kind + 'NativeOwnedShader'] = pin(R / modules[kind]['module'])
header = (C / 'retained_kernels.h').read_text()
available = {stage: re.search(r'k' + prefix + r'FmaAvailable = (true|false);', header)[1] == 'true'
             for stage, prefix in [('transport', 'Transport'), ('endpoint', 'Endpoint')]}
creation, = doc(W.parent / 'provider-history/candidate-default/device-features.json')['device_creations']
props = doc(W.parent / 'provider-history/candidate-default/float-controls-properties.json')
assert props and all(x['values']['shaderSignedZeroInfNanPreserveFloat32'] == props[0]['values']['shaderSignedZeroInfNanPreserveFloat32'] for x in props)
identity = dict(control8['record']['identity'])
identity['fma_fp32_enabled'] = str(int(any(x['fma16_32_64'][1] == 1 for x in creation['feature_chain'])))
identity['preserves_fp32_signed_zero_inf_nan'] = str(props[0]['values']['shaderSignedZeroInfNanPreserveFloat32'])
context = SimpleNamespace(value={'native_identity': identity},
    arms={'candidate': {'arrays': arrays, 'fma_available': available}})
original_select = old.selected_names
old.selected_names = lambda context, arm, wide=False: original_select(context, arm, wide) if wide else native
accepted = {}
for i in range(1, 8):
    label = 'control-' + str(i); folder = RUN / label
    owner = doc(folder / 'owner.json'); action = actions[label]
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert owner['source'] == {'revision': action['revision'], 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    assert owner['wall_seconds'] <= 300 and owner['peak_sampled_rss_kib'] <= 4194304
    assert owner['minimum_sampled_available_kib'] >= 2097152
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
    before = doc(folder / 'inputs-before.json'); assert before == doc(folder / 'inputs-after.json')
    for key, value in before.items():
        p = Path(key) if key.startswith('/') else R / key
        if key == '/usr/lib/x86_64-linux-gnu/liblzma.so.5.4.5':
            # The OS upgraded this one file after every owned run closed.
            # Require the exact previously frozen bytes, privately recovered;
            # preserve the original before/after seals and old path verbatim.
            restoration = doc(W.parent / 'library-restoration/identity-restoration.json')
            assert restoration['historical_path'] == key
            assert restoration['historical_identity'] == {k: value[k] for k in ['bytes', 'sha256']}
            p = R / restoration['private_restoration']
        assert pin(p) == {k: value[k] for k in ['bytes', 'sha256']}, key
    env = owner['declared_output_environment']
    assert env == {'GTEST_OUTPUT': 'xml:' + str(folder / 'gtest.xml'),
        'VK_DRIVER_FILES': provider['ICD'], 'SIRIUS_VULKAN_DEVICE': '0', 'MESA_SHADER_CACHE_DISABLE': 'true'}
    assert owner['provider_mappings'] and all(x['owned_session'] == owner['child_pid']
        and x['maps'] and all(provider['library'] in line for line in x['maps']) for x in owner['provider_mappings'])
    props = old.xml_properties(folder / 'gtest.xml', label)
    row = {'passed': True, 'case': action['case'], 'owner': pin(folder / 'owner.json'),
        'XML': pin(folder / 'gtest.xml'), 'wall_seconds': owner['wall_seconds'],
        'peak_sampled_rss_kib': owner['peak_sampled_rss_kib'], 'properties': props}
    if i == 1:
        row['complete_native_routes'] = old.factored_routes(
            {'props': props, 'folder': folder, 'arm': 'candidate'}, context)
        assert row['complete_native_routes']['routes_count'] == 2
        assert row['complete_native_routes']['complete_readbacks'] == 28
    accepted[label] = row
out = {'all8_original_controls_pass': True, 'all22_original_benefit_checks_pass': True,
       'control8_acceptance': pin(W / 'control-8-acceptance.json'),
       'benefit': pin(W / 'benefit.json'), 'decimal': pin(W / 'decimal-agreement.json'),
       'remaining7': accepted, 'reader': pin(Path(__file__).resolve()),
       'postexecution_library_identity_restoration': doc(W.parent / 'library-restoration/identity-restoration.json'),
       'scope': 'Original eight finite controls and sole benefit gate only. Does not admit a cold frame, full scientific/Mandatory/platform/native-window/installed-release qualification.'}
p = W / 'controls-1-8-acceptance.json'; assert not p.exists()
p.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({'all8_original_controls_pass': True, 'all22_original_benefit_checks_pass': True,
                  'factored_actual_precision_routes': 2, 'factored_complete_readbacks': 28}))
