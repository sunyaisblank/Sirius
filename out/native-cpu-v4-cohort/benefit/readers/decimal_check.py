from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-cohort/benefit/readers/decimal_check.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='c1861053cb6216ef91f13852b3f3447d03a8d665ba30f3f0c1bf803502bd6ac3'
source=basis.read_text()
assert source.count('out/native-cpu-local-identities-cohort')==1
source=source.replace('out/native-cpu-local-identities-cohort','out/native-cpu-v4-cohort')
exec(compile(source,str(basis),"exec"))
