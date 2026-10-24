"""Reuse the existing incremental build planner and resource limits."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-slang-bridge/prepare_provider_build.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '5a91184fb45804da505f76552448fdb617a8eff6220909c5351bb031eee531d1'
exec(compile(basis.read_text(), str(basis), 'exec'))
