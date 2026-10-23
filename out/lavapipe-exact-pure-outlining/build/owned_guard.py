from pathlib import Path
P = Path.cwd() / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py'
exec(compile(P.read_text(), str(P), 'exec'))
