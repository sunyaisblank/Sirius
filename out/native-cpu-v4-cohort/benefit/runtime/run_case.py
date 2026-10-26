from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-cohort/benefit/runtime/run_case.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='286d5cfa6c9644e374b147595d4ad4d2b31d9c3b68e0350749dbb1fcf7364912'
source=basis.read_text()
assert source.count('out/native-cpu-local-identities-cohort')==1
source=source.replace('out/native-cpu-local-identities-cohort','out/native-cpu-v4-cohort')
exec(compile(source,str(basis),"exec"))
