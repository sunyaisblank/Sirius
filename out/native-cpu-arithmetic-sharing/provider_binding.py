"""Bind original controls to actual current modules, provider and matched host."""
from pathlib import Path
import hashlib
import json

R = Path.cwd()
W = Path(__file__).resolve().parent
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
PRIOR = R / 'out/native-original-input-ownership'
def doc(p): return json.loads(p.read_text())
def pin(p):
    b = p.read_bytes()
    return {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
def terminal(p):
    o = doc(p)
    assert o['passed'] and o['status'] == 'terminal' and o['child_exit'] == 0
    assert o['source_unchanged'] and o['whole_inputs_unchanged']
    assert not o['remaining_owned_processes'] and not o['cleanup_errors']
    assert o['observed_births'] and all(b['absent'] for b in o['observed_births'])
    assert doc(p.parent / 'inputs-before.json') == doc(p.parent / 'inputs-after.json')
    return o

def admit():
    review = doc(W / 'provider/review-admission.json')
    assert review['passed']
    assert pin(W / 'provider/provider.patch') == review['prototype_patch']
    assert pin(W / 'provider/source-contract.json') == review['source_contract']
    for name, expected in review['source_files'].items():
        assert pin(S / name) == expected
    assert doc(W / 'emission-result.json')['passed']
    terminal(W / 'build/modules/owner.json')
    terminal(W / 'build/matched/owner.json')
    terminal(W / 'provider/build/owner.json')
    matched = doc(W / 'matched-result.json')
    assert matched['passed'] and matched['same_actual_executable_for_both_provider_arms']
    assert matched['actual_three_module_spans_present_once_in_selector_object_and_final_ELF']
    exe = R / matched['executable']
    assert pin(exe)['sha256'] == matched['executable_sha256']
    original_selector = doc(W / 'source-scope.json')['host_selector']
    assert pin(W / 'retained_compute.cpp')['sha256'] == original_selector['sha256']
    providers = doc(W / 'providers.json')['providers']
    provider = providers['candidate']
    for p in providers.values():
        assert pin(Path(p['library'])) == p['library_identity']
        assert pin(Path(p['ICD'])) == p['ICD_identity']
    built = doc(W / 'provider/provider-build-result.json')
    assert built['passed'] and built['library_identity'] == provider['library_identity']
    # Primitive proof supplies unchanged source context, never changed-module numerical results.
    for label in ['candidate-default', 'candidate-poison']:
        proof = doc(PRIOR / 'reused/prior-proof/lavapipe-ieee32-validation' / label / 'numerical-and-boundary-result.json')
        assert proof['numerical_pass'] and proof['all_observed_callback_states_restored']
    expected_props = doc(PRIOR / 'provider-history/candidate-default/float-controls-properties.json')
    expected_features, = doc(PRIOR / 'provider-history/candidate-default/device-features.json')['device_creations']
    modules = doc(W / 'compiler-result.json')['modules']
    code = doc(W / 'ordinary/ordinary-code-review.json')
    assert code['structural_premise_observed']
    binary = exe.read_bytes()
    for row in modules.values():
        module = R / row['module']
        assert pin(module) == {k: row[k] for k in ['bytes', 'sha256']}
        assert binary.count(module.read_bytes()) == 1
    for label, kind in [('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
        folder = W / 'ordinary' / label
        o = terminal(folder / 'owner.json')
        v = doc(folder / 'capture-verification.json')
        assert v['passed'] and v['actual_provider_sha256'] == provider['library_identity']['sha256']
        assert v['actual_complete_module']['sha256'] == modules[kind]['sha256']
        assert v['callback_restoration_observed'] and v['actual_entry']['actual_mxcsr'] & 0xe040 == 0
        row = code['rows'][kind]
        assert row['structural_premise_observed'] and row['complete_bulk_clear_observed']
        assert row['original_Default_machine_pipeline_observed']
        assert row['total_called_function_bytes'] < row['prior_best_called_text_bytes']
        assert o['provider_mappings'] and all(provider['library'] in line
            for m in o['provider_mappings'] for line in m['maps'])
        props = doc(folder / 'float-controls-properties.json')
        assert props and all(p['values'] == expected_props[0]['values'] for p in props)
        features, = doc(folder / 'device-features.json')['device_creations']
        assert features['features']['values'] == expected_features['features']['values']
        assert [x['fma16_32_64'] for x in features['feature_chain']] == [x['fma16_32_64'] for x in expected_features['feature_chain']]
    pipeline = (S / 'src/gallium/frontends/lavapipe/lvp_pipeline.c').read_text()
    assert '0x27,0xc1,0xce,0xb7' not in pipeline and '0xc0,0xf3,0x8b,0x8b' not in pipeline
    assert modules['Transport']['sha256'] == '36b5e823b363f68c519cb10073cc113899d0183be3314434b003341d3531088f'
    return True
