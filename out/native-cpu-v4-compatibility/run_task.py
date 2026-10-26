from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-compatibility/run_task.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='e63294a00c9d79a7403a96d52415feeff39739b41d7e772facc74575f2028585'
source=basis.read_text()
assert source.count('exact OLD/new S raw Transport code')==1
source=source.replace('exact OLD/new S raw Transport code','exact original-Q/unchanged-AN raw Transport with new guarded optional math composition',1)
exec(compile(source,str(basis),"exec"))
