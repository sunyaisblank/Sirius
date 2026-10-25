from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-row-ownership/cold-frame/run_owned.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='34f021ffe2064fe489181e2a370f33caf8c3a75ee60dd5088629fee9a68a8aab'
exec(compile(basis.read_text(),str(basis),'exec'))
