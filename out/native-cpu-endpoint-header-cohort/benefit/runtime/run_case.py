from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from authority import load
exec(compile(load('runtime/run_case.py'), '/home/astra/.project/Sirius/out/native-cpu-v4-cohort/benefit/runtime/run_case.py', "exec"))
