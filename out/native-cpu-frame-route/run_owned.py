from pathlib import Path
basis=Path.cwd()/'out/native-cpu-row-ownership/cold-frame/run_owned.py'
exec(compile(basis.read_text(),str(basis),'exec'))
