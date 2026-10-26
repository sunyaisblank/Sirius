from pathlib import Path
import hashlib
basis=Path.cwd()/'out/native-cpu-slang-bridge/prepare_provider_build.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='5a91184fb45804da505f76552448fdb617a8eff6220909c5351bb031eee531d1'
source=basis.read_text()
exec(compile(source,str(basis),"exec"))
