from pathlib import Path
basis=Path.cwd()/'out/native-cpu-row-ownership/cold-frame/owned_guard.py'
exec(compile(basis.read_text(),str(basis),'exec'))
