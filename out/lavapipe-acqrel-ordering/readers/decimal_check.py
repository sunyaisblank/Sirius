from pathlib import Path
import hashlib
_basis = Path.cwd() / 'out/native-original-input-ownership/readers/decimal_check.py'
assert hashlib.sha256(_basis.read_bytes()).hexdigest() == 'd27a3e94775edbdba42b48569d0e09a794b720fc13bc94d92787f9c62b1f305b'
_source = _basis.read_text()
assert _source.count("RUN=R/'out/native-original-input-ownership/runtime'") == 1
_source = _source.replace("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/lavapipe-acqrel-ordering/runtime'")
exec(compile(_source, str(_basis), "exec"))
