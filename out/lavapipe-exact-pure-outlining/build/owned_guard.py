from pathlib import Path
P = Path.cwd() / 'out/lavapipe-acqrel-ordering/runtime/owned_guard.py'
exec(compile(P.read_text(), str(P), 'exec'))
