from pathlib import Path
p=Path.cwd()/'out/native-cpu-row-ownership/cold-frame/owned_guard.py'
exec(compile(p.read_text(),str(p),"exec"))
