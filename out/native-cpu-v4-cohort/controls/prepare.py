from pathlib import Path
from control_authority import load
source=load('prepare.py')
anchor="assert doc(current / 'control-8-acceptance.json')['passed']"
assert source.count(anchor)==1
source=source.replace(anchor,anchor+"\nassert doc(current / 'control-8-acceptance.json')['optional_math']['all_selected_optional_math_observed']",1)
exec(compile(source,'/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/controls/prepare.py',"exec"))
