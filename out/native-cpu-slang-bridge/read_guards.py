"""Compose actual native refusal, unchanged canary, ordinary JIT and FP envelopes."""
from pathlib import Path
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent
O = W / 'ordinary'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


actions = doc(O / 'actions.json')
build = doc(W / 'provider/provider-build-result.json')
linked = doc(W / 'provider/linked-native-proof.json')
assert linked['passed'] and linked['library_identity'] == build['library_identity']
assert pin(R / build['library']) == build['library_identity']
modules = doc(R / 'out/native-cpu-row-ownership/compiler-result.json')['modules']
prior = R / 'out/native-original-input-ownership/provider-history/candidate-default'
expected_properties = doc(prior / 'float-controls-properties.json')[0]['values']
expected_features, = doc(prior / 'device-features.json')['device_creations']
rows = {}
for slug, stage in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
    folder = O / slug
    owner = doc(folder / 'owner.json')
    assert owner['passed'] and owner['status'] == 'terminal' and owner['child_exit'] == 0
    assert owner['source'] == {'revision': actions[slug]['revision'], 'status': ''}
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors'] and owner['stop_reason'] is None
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    before = doc(folder / 'inputs-before.json')
    assert before == doc(folder / 'inputs-after.json')
    for name, identity in before.items():
        assert pin(R / name) == {key: identity[key] for key in ('bytes', 'sha256')}, name
    assert owner['provider_mappings'] and all(build['library'] in line for mapping in owner['provider_mappings'] for line in mapping['maps'])
    native, = doc(folder / 'native-calls.json')
    verification = doc(folder / 'capture-verification.json')
    assert native['function'] == 'sirius_native_' + slug + '_try_packet'
    assert native['library'] == str(R / build['library'])
    assert native['returned'] and native['handled'] is False
    assert native['input_bytes'] == native['output_bytes'] == 4
    assert native['output_before'] == native['output_after'] == 'baaadbda'
    assert native['first_x'] == native['y'] == native['z'] == 0
    assert native['grid_x'] == native['grid_y'] == native['grid_z'] == native['groups'] == 1
    assert native['mxcsr'] & 0xe040 == 0
    assert verification['passed'] and verification['callback_restoration_observed']
    assert verification['actual_provider_sha256'] == build['library_identity']['sha256']
    assert verification['actual_complete_module']['sha256'] == modules[stage]['sha256']
    assert verification['actual_complete_module']['actual_span_matches_frozen_input']
    entry = verification['actual_entry']
    assert entry['worker_thread'] == native['thread']
    assert entry['actual_mxcsr'] == native['mxcsr']
    boundary, = doc(folder / 'compute-boundaries.json')
    assert boundary['returned'] and boundary['worker'] == native['thread']
    assert boundary['mxcsr_before'] == boundary['mxcsr_after']
    properties = doc(folder / 'float-controls-properties.json')
    assert properties and all(row['values'] == expected_properties for row in properties)
    features, = doc(folder / 'device-features.json')['device_creations']
    assert features['features']['values'] == expected_features['features']['values']
    assert [row['fma16_32_64'] for row in features['feature_chain']] == [row['fma16_32_64'] for row in expected_features['feature_chain']]
    rows[stage] = {'passed': True, 'owner': pin(folder / 'owner.json'),
                   'native_call': pin(folder / 'native-calls.json'), 'native_returned_false': True,
                   'output_unchanged': True, 'same_worker_real_JIT_called': True,
                   'callback_full_MXCSR_restored': True, 'features_properties_match_prior': True,
                   'wall_seconds': owner['wall_seconds'], 'peak_sampled_rss_kib': owner['peak_sampled_rss_kib']}
record = {'passed': True, 'rows': rows, 'provider_identity': build['library_identity'],
          'static_native_proof': pin(W / 'provider/linked-native-proof.json'),
          'scope': 'Three finite original4B negative boundary captures only: actual linked native entry returnsfalse beforewrites, sameworker executes ordinary JIT, completecanary holds and FPstate restores. Does not exercise successful native arithmetic. Fresh originalC8 must prove successful native admission/all32 exactoutputs before soleABBA22+Decimal/all8/frame. No performance/numerical/scientific/release adoption.'}
assert not (O / 'guard-acceptance.json').exists()
(O / 'guard-acceptance.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'negative_guard_stages_pass': len(rows), 'native_arithmetic_proved': False}))
