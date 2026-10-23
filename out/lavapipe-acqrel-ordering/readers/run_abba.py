from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-original-input-ownership/readers/run_abba.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == '14177971258ceb4a269549f1c0f28055dec6d165f83e46f6047953e53be0d7eb'
_source = _basis.read_text()
assert _source.count("RUN=R/'out/native-original-input-ownership/runtime'") == 1
_source = _source.replace("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/lavapipe-acqrel-ordering/runtime'")
assert _source.count("RUN/'run_case.py'") == 1
_source = _source.replace("RUN/'run_case.py'", "RUN.parent/'run_case.py'")
exec(compile(_source, str(_basis), "exec"))
