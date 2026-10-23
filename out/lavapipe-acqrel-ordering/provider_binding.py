"""Compose original precision proof with actual new provider code and observations."""
from pathlib import Path
import hashlib, json

R = Path.cwd(); W = R / 'out/lavapipe-acqrel-ordering'
OLD = R / 'out/native-original-input-ownership'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def admit():
    matched = doc(OLD / 'matched-result.json')
    assert matched['passed'] and matched['same_actual_executable_for_both_provider_arms']
    assert pin(R / matched['executable'])['sha256'] == matched['executable_sha256']
    source = doc(W / 'source-admission.json'); assert source['passed']
    assert len(source['changed_files']) == 6
    provider = doc(W / 'providers.json')['providers']['candidate']
    assert pin(Path(provider['library'])) == provider['library_identity']
    assert pin(Path(provider['ICD'])) == provider['ICD_identity']
    build = doc(W / 'build/owner.json')
    assert build['passed'] and build['source_unchanged'] and build['whole_inputs_unchanged']
    assert build['status'] == 'terminal' and not build['remaining_owned_processes'] and not build['cleanup_errors']
    for name in ['candidate-default', 'candidate-poison']:
        proof = doc(OLD / 'reused/prior-proof/lavapipe-ieee32-validation' / name / 'numerical-and-boundary-result.json')
        assert proof['numerical_pass'] and proof['all_observed_callback_states_restored']
    old_properties = doc(OLD / 'provider-history/candidate-default/float-controls-properties.json')
    old_creation, = doc(OLD / 'provider-history/candidate-default/device-features.json')['device_creations']
    code = doc(W / 'ordinary/ordinary-code-review.json'); assert code['structural_premise_observed']
    compiler = doc(OLD / 'compiler-result.json')
    for label, kind in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
        folder = W / 'ordinary' / label
        owner = doc(folder / 'owner.json')
        assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        assert owner['status'] == 'terminal' and not owner['remaining_owned_processes'] and not owner['cleanup_errors']
        actual = doc(folder / 'capture-verification.json')
        assert actual['passed'] and actual['actual_provider_sha256'] == provider['library_identity']['sha256']
        assert actual['actual_complete_module']['sha256'] == compiler['modules'][kind]['sha256']
        assert actual['callback_restoration_observed'] and actual['actual_entry']['actual_mxcsr'] & 0xe040 == 0
        assert owner['provider_mappings'] and all(provider['library'] in line
            for mapping in owner['provider_mappings'] for line in mapping['maps'])
        properties = doc(folder / 'float-controls-properties.json')
        assert properties and all(x['values'] == old_properties[0]['values'] for x in properties)
        creation, = doc(folder / 'device-features.json')['device_creations']
        assert creation['features']['values'] == old_creation['features']['values']
        assert [x['fma16_32_64'] for x in creation['feature_chain']] == [x['fma16_32_64'] for x in old_creation['feature_chain']]
        assert code['rows'][kind]['structural_premise_observed']
    return True
