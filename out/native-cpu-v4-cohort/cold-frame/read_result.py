"""Original frame/XML/owner/resource reader plus current optional-C8 binding."""
from pathlib import Path
from frame_authority import load
source=load('read_result.py')
anchor='source = basis.read_text()'
assert source.count(anchor)==1
injection='source = basis.read_text()\nanchor="assert doc(current / \'control-8-acceptance.json\')[\'passed\']"\nassert source.count(anchor)==1\nsource=source.replace(anchor,anchor+"\\nassert doc(current / \'control-8-acceptance.json\')[\'optional_math\'][\'all_selected_optional_math_observed\']",1)\n'
source=source.replace(anchor,injection,1)
exec(compile(source,'/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/cold-frame/read_result.py',"exec"))
