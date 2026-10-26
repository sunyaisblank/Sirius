from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-normalization-tags/apply_provider.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='21b6f4f7f86f54b65407907ce63c6e096529ae89e3e084941b1cbdb8cc7cfda2'
source=basis.read_text()
assert source.count('AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback')==1
source=source.replace('AJ caches the exact private-register normalization predicate with complete write invalidation and unchanged arithmetic/fallback','AO composes separately guarded same-eight-lane optional v4 math and byte-preserved v3 fallback; original eligibility/FP scopes/raw3/physical contracts remain',1)
exec(compile(source,str(basis),"exec"))
