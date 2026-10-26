from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-compatibility/link.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='f5d43a16506d29ba88ebac641702c1cf580c6e434feb38ee2d9a7e0c3c8ec548'
source=basis.read_text()
exec(compile(source,str(basis),"exec"))
