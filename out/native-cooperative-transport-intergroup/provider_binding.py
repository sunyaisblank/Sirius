"""Compose unchanged ordinary Transport proof with exact current E/DP code proof."""
from pathlib import Path
import hashlib, json

R = Path.cwd(); W = Path(__file__).resolve().parent
CODE = R / 'out/lavapipe-intergroup-row-owner'
OLD = R / 'out/native-original-input-ownership'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def terminal(p):
    owner = doc(p)
    assert owner['passed'] and owner['status'] == 'terminal'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
    return owner
def admit():
    source = doc(W / 'source-scope.json'); assert source['source_review_passed']
    assert pin(W / 'retained_compute.cpp') == source['host_selector']
    assert pin(W / 'headers/retained_kernels.h') == source['header']
    terminal(W / 'matched/owner.json')
    matched = doc(W / 'matched-result.json')
    assert matched['passed'] and matched['same_actual_executable_for_both_provider_arms']
    assert matched['actual_three_module_spans_present_once_in_selector_object_and_final_ELF']
    assert pin(R / matched['executable'])['sha256'] == matched['executable_sha256']
    provider = doc(W / 'providers.json')['providers']['candidate']
    assert pin(Path(provider['library'])) == provider['library_identity']
    assert pin(Path(provider['ICD'])) == provider['ICD_identity']
    assert provider['library_identity']['sha256'] == source['provider_sha256']
    terminal(CODE / 'build/owner.json')
    built = doc(CODE / 'provider-build-result.json')
    assert built['passed'] and built['library_identity'] == provider['library_identity']
    reuse = doc(W / 'code-reuse.json'); assert reuse['source_composition_admitted']
    for key in ['Transport', 'Endpoint_DopriPhase']:
        proof = reuse[key]; assert pin(R / proof['proof']) == {k: proof[k] for k in ['bytes', 'sha256']}
    compiler = doc(W / 'compiler-result.json')['modules']
    old_code = doc(R / reuse['Transport']['proof'])
    assert old_code['structural_premise_observed']
    assert old_code['rows']['Transport']['complete_module']['sha256'] == compiler['Transport']['sha256']
    # This is reuse of an unchanged source/lowering path, not a claim that
    # Transport was recaptured under this provider. Fresh C8 must bind its
    # actual selected module, residency and complete readbacks.
    assert hashlib.sha1((R / compiler['Transport']['module']).read_bytes()).hexdigest() == reuse['Transport']['whitelist_sha1_absent']
    for label in ['candidate-default', 'candidate-poison']:
        proof = doc(OLD / 'reused/prior-proof/lavapipe-ieee32-validation' / label / 'numerical-and-boundary-result.json')
        assert proof['numerical_pass'] and proof['all_observed_callback_states_restored']
    properties = doc(OLD / 'provider-history/candidate-default/float-controls-properties.json')
    creation, = doc(OLD / 'provider-history/candidate-default/device-features.json')['device_creations']
    code = doc(R / reuse['Endpoint_DopriPhase']['proof']); assert code['structural_premise_observed']
    for label, kind in [('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
        folder = CODE / 'ordinary' / label; owner = terminal(folder / 'owner.json')
        actual = doc(folder / 'capture-verification.json')
        assert actual['passed'] and actual['actual_provider_sha256'] == provider['library_identity']['sha256']
        assert actual['actual_complete_module']['sha256'] == compiler[kind]['sha256']
        assert code['rows'][kind]['structural_premise_observed']
        assert actual['callback_restoration_observed'] and actual['actual_entry']['actual_mxcsr'] & 0xe040 == 0
        assert owner['provider_mappings'] and all(provider['library'] in line
            for mapping in owner['provider_mappings'] for line in mapping['maps'])
        observed = doc(folder / 'float-controls-properties.json')
        assert observed and all(x['values'] == properties[0]['values'] for x in observed)
        current, = doc(folder / 'device-features.json')['device_creations']
        assert current['features']['values'] == creation['features']['values']
        assert [x['fma16_32_64'] for x in current['feature_chain']] == [x['fma16_32_64'] for x in creation['feature_chain']]
    return True
