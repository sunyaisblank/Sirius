from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-cohort/benefit/readers/native_evidence.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='bf50b12e63d42c0715007e173e25ea2abbf01f71030b9d62099be034d59e40e8'
source=basis.read_text()
exec(compile(source,str(basis),"exec"))
