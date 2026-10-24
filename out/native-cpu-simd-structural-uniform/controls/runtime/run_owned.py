from pathlib import Path
import hashlib
source = Path.cwd() / 'out/native-cpu-row-ownership/cold-frame/run_owned.py'
assert hashlib.sha256(source.read_bytes()).hexdigest() == '34f021ffe2064fe489181e2a370f33caf8c3a75ee60dd5088629fee9a68a8aab'
exec(compile(source.read_text(), str(source), 'exec'))
