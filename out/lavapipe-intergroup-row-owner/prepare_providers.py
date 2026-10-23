"""Bind actual staged provider bytes; historical records supply expectations only."""
from pathlib import Path
import hashlib, json

R = Path.cwd(); W = Path(__file__).resolve().parent
def doc(p): return json.loads(p.read_text())
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
assert doc(W / 'ordinary/ordinary-code-review.json')['structural_premise_observed']
result = doc(W / 'provider-build-result.json'); assert result['passed']
old = doc(R / 'out/native-original-input-ownership/providers.json')['providers']
baseline = dict(old['baseline'])
assert pin(Path(baseline['library'])) == baseline['library_identity']
assert doc(Path(baseline['ICD'])) == {'file_format_version': '1.0.0',
    'ICD': {'library_path': baseline['library'], 'api_version': '1.4.305'}}
baseline['historical_ICD_identity_expected_only'] = baseline['ICD_identity']
baseline['ICD_identity'] = pin(Path(baseline['ICD']))
baseline['descriptor_scope'] = 'Seal the actual existing descriptor pointing at the exact original stock library. Historical descriptor digest differs; no descriptor or provider byte rewritten.'
candidate = {'library': str(R / result['library']), 'ICD': str(R / result['ICD']),
    'library_identity': result['library_identity'], 'ICD_identity': result['ICD_identity'],
    'historical_library': old['candidate']['historical_library'],
    'historical_identity_only': 'Historical hardware/precision expectations only. Actual new provider bytes, all-three current code/features/properties/MXCSR and fresh original runtime identities/readbacks are separately required.'}
assert pin(Path(candidate['library'])) == candidate['library_identity']
assert pin(Path(candidate['ICD'])) == candidate['ICD_identity']
target = W / 'providers.json'; assert not target.exists()
target.write_text(json.dumps({'providers': {'baseline': baseline, 'candidate': candidate}}, indent=2) + '\n')
print(json.dumps({'exact_stock_baseline_reused': True, 'actual_staged_provider': candidate['library_identity'],
                  'new_CPP_probe_shader_compilations': 0, 'runtime_admitted': False}))
