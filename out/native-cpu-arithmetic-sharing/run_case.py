from pathlib import Path
import hashlib
p=Path.cwd()/'out/lavapipe-acqrel-ordering/run_case.py'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='f9c87d6108ffd8ce17d2cf9867e0ab6340baef65310b6437eb6407a321c5b966'
s=p.read_text()
old="bin/linux-gcc/native-original-input-ownership/matched/audited_backend_tests"
assert s.count(old)==1
s=s.replace(old,"bin/linux-gcc/native-cpu-arithmetic-sharing/matched/audited_backend_tests")
exec(compile(s,str(p),"exec"))
