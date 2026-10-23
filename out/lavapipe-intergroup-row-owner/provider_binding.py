"""Reuse the original acceptance authority with exact current artifact paths."""
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/lavapipe-acqrel-ordering/provider_binding.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'cc18e9b6dc621d9bb9d34ed6ec163a0af2f39f99a00e9ad175e0a7711dc68b51'
_source = _basis.read_text()
assert 'out/lavapipe-acqrel-ordering' in _source
_source = _source.replace('out/lavapipe-acqrel-ordering', 'out/lavapipe-intergroup-row-owner')
assert "matched = doc(OLD / 'matched-result.json')" in _source
_source = _source.replace("matched = doc(OLD / 'matched-result.json')", "matched = doc(R / 'out/native-cpu-row-ownership/matched-result.json')")
assert "doc(W / 'source-admission.json')" in _source
_source = _source.replace("doc(W / 'source-admission.json')", "doc(W / 'review-admission.json')")
assert "assert len(source['changed_files']) == 6" in _source
_source = _source.replace("assert len(source['changed_files']) == 6", "assert source['reverse_patch_recovers_exact_prior_six_file_basis'] and len(source['changed_from_prior_trial']) == 7")
assert "compiler = doc(OLD / 'compiler-result.json')" in _source
_source = _source.replace("compiler = doc(OLD / 'compiler-result.json')", "compiler = doc(R / 'out/native-cpu-row-ownership/compiler-result.json')")
exec(compile(_source, str(_basis), "exec"))
