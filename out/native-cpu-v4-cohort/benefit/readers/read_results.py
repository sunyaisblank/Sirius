from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-local-identities-cohort/benefit/readers/read_results.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='3f3f1a7e3703a084460ad077777ebfd4defa5a1290fe9586ef5f664765dfe07b'
source=basis.read_text()
assert source.count('out/native-cpu-local-identities-cohort')==1
source=source.replace('out/native-cpu-local-identities-cohort','out/native-cpu-v4-cohort')
exec(compile(source,str(basis),"exec"))
