from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from control_authority import load
source=load('runtime/run_owned.py')
exec(compile(source,'/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/controls/runtime/run_owned.py',"exec"))
