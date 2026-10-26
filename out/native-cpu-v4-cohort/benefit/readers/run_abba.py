from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-cohort/benefit/readers/run_abba.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='a5b93d12093307cb8f80150c6b74ad2364455f4a697e871417f49ee695ba90a9'
source=basis.read_text()
assert source.count('out/native-cpu-local-identities-cohort')==1
source=source.replace('out/native-cpu-local-identities-cohort','out/native-cpu-v4-cohort')
exec(compile(source,str(basis),"exec"))
