"""Admit actual provider/code identities before invoking any original workload."""
from pathlib import Path
import hashlib, json
R = Path.cwd(); W = Path(__file__).resolve().parent
OLD = R / 'out/native-cooperative-transport-intergroup'
PRIOR = R / 'out/native-original-input-ownership'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size, 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
def terminal(p):
    owner = doc(p)
    assert owner['passed'] and owner['status'] == 'terminal'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert owner['observed_births'] and all(x['absent'] for x in owner['observed_births'])
    return owner

def admit():
    review = doc(W / 'review-admission.json'); assert review['passed']
    assert pin(W / 'machine-codegen.patch') == review['machine_codegen_patch']
    assert pin(W / 'source-contract.json') == review['source_contract']
    for name, expected in review['source_files'].items(): assert pin(S / name) == expected, name
    source = doc(OLD / 'source-scope.json'); assert source['source_review_passed']
    assert pin(OLD / 'retained_compute.cpp') == source['host_selector']
    assert pin(OLD / 'headers/retained_kernels.h') == source['header']
    terminal(OLD / 'matched/owner.json')
    matched = doc(W / 'matched-result.json')
    assert matched['passed'] and matched['same_actual_executable_for_both_provider_arms']
    assert matched['actual_three_module_spans_present_once_in_selector_object_and_final_ELF']
    assert pin(R / matched['executable'])['sha256'] == matched['executable_sha256'] == 'f958a60563b64c79d3f2fb6da2815bde0effdbbe49f121034cfd212ab17b27ff'
    provider = doc(W / 'providers.json')['providers']['candidate']
    assert pin(Path(provider['library'])) == provider['library_identity']
    assert pin(Path(provider['ICD'])) == provider['ICD_identity']
    terminal(W / 'build/owner.json')
    built = doc(W / 'provider-build-result.json')
    assert built['passed'] and built['library_identity'] == provider['library_identity']
    for label in ['candidate-default', 'candidate-poison']:
        proof = doc(PRIOR / 'reused/prior-proof/lavapipe-ieee32-validation' / label / 'numerical-and-boundary-result.json')
        assert proof['numerical_pass'] and proof['all_observed_callback_states_restored']
    properties = doc(PRIOR / 'provider-history/candidate-default/float-controls-properties.json')
    creation, = doc(PRIOR / 'provider-history/candidate-default/device-features.json')['device_creations']
    compiler = doc(W / 'compiler-result.json')['modules']
    code = doc(W / 'ordinary/ordinary-code-review.json')
    assert code['structural_premise_observed']
    for label, kind in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
        module = R / compiler[kind]['module']
        assert pin(module) == {k: compiler[kind][k] for k in ['bytes', 'sha256']}
        assert (R / matched['executable']).read_bytes().count(module.read_bytes()) == 1
        folder = W / 'ordinary' / label; owner = terminal(folder / 'owner.json')
        actual = doc(folder / 'capture-verification.json')
        assert actual['passed'] and actual['actual_provider_sha256'] == provider['library_identity']['sha256']
        assert actual['actual_complete_module']['sha256'] == compiler[kind]['sha256']
        assert code['rows'][kind]['structural_premise_observed']
        assert code['rows'][kind]['machine_policy_isolation_observed']
        assert code['rows'][kind]['optimized_IR_unchanged']
        assert actual['callback_restoration_observed'] and actual['actual_entry']['actual_mxcsr'] & 0xe040 == 0
        assert owner['provider_mappings'] and all(provider['library'] in line
            for mapping in owner['provider_mappings'] for line in mapping['maps'])
        observed = doc(folder / 'float-controls-properties.json')
        assert observed and all(x['values'] == properties[0]['values'] for x in observed)
        current, = doc(folder / 'device-features.json')['device_creations']
        assert current['features']['values'] == creation['features']['values']
        assert [x['fma16_32_64'] for x in current['feature_chain']] == [x['fma16_32_64'] for x in creation['feature_chain']]
    return True
