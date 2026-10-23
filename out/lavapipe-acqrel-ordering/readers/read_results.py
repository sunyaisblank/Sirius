from pathlib import Path
import importlib.util
_binding = Path(__file__).resolve().parents[1] / 'provider_binding.py'
_spec = importlib.util.spec_from_file_location('actual_provider_binding', _binding)
_module = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_module)
_module.admit()
from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-original-input-ownership/readers/read_results.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'e7ef5a184fd2f25f5d0d606baadcae0e8ec8fb2f4e14e88c41e24915b556034e'
_source = _basis.read_text()
assert _source.count("RUN=R/'out/native-original-input-ownership/runtime'") == 1
_source = _source.replace("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/lavapipe-acqrel-ordering/runtime'")
assert _source.count("doc(W.parent/'compiler-result.json')") == 1
_source = _source.replace("doc(W.parent/'compiler-result.json')", "doc(R/'out/native-original-input-ownership/compiler-result.json')")
assert _source.count("doc(W.parent/'previous-control-8-acceptance.json')") == 1
_source = _source.replace("doc(W.parent/'previous-control-8-acceptance.json')", "doc(R/'out/native-original-input-ownership/readers/control-8-acceptance.json')")
exec(compile(_source, str(_basis), "exec"))
