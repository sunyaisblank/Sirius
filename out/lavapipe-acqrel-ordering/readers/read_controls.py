from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-original-input-ownership/readers/read_controls.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '2c645bf4033da6a31cb0bee68cdf1e1d44fba07d371766168b1646f707e15f65'
_source = _basis.read_text()
assert _source.count("doc(W.parent / 'compiler-result.json')") == 1
_source = _source.replace("doc(W.parent / 'compiler-result.json')", "doc(R / 'out/native-original-input-ownership/compiler-result.json')")
assert _source.count("doc(W.parent / 'provider-history/candidate-default/device-features.json')") == 1
_source = _source.replace("doc(W.parent / 'provider-history/candidate-default/device-features.json')", "doc(W.parent / 'ordinary/transport/device-features.json')")
assert _source.count("doc(W.parent / 'provider-history/candidate-default/float-controls-properties.json')") == 1
_source = _source.replace("doc(W.parent / 'provider-history/candidate-default/float-controls-properties.json')", "doc(W.parent / 'ordinary/transport/float-controls-properties.json')")
exec(compile(_source, str(_basis), "exec"))
