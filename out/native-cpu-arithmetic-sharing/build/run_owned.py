from pathlib import Path
p=Path.cwd()/'out/native-cpu-row-ownership/cold-frame/run_owned.py'
exec(compile(p.read_text(),str(p),"exec"))
