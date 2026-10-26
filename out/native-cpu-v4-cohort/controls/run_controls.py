from pathlib import Path
from control_authority import load
source=load('run_controls.py')
exec(compile(source,'/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/controls/run_controls.py',"exec"))
