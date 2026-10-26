from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-cohort/benefit/readers/numeric_gates.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='5317872eae1fccf0ba9a04f552f3c303e309fda2f0279e573ef1137d8772ca12'
source=basis.read_text()
exec(compile(source,str(basis),"exec"))
