from pathlib import Path
import hashlib
basis = Path.cwd() / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py'
data = basis.read_bytes()
assert hashlib.sha256(data).hexdigest() == '747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056'
exec(compile(data, str(basis), 'exec'))
