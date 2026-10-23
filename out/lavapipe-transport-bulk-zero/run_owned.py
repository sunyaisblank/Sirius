from pathlib import Path
P=Path(__file__).resolve().parents[1]/'native-cpu-row-ownership/cold-frame/run_owned.py'
exec(compile(P.read_text(),str(P),'exec'))
