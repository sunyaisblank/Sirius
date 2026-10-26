from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from control_authority import load
source=load('readers/read_controls.py')
anchor="exec(compile(source, str(basis), 'exec'))"
assert source.count(anchor)==1
addition='anchor="p = W / \'controls-1-8-acceptance.json\'"\nassert source.count(anchor)==1\nextra="assert control8[\'optional_math\'][\'all_selected_optional_math_observed\']\\nout[\'optional_C8_runtime_coverage\'] = control8[\'optional_math\']\\nout[\'provider_scope\'] = \'Same original H host/AO guarded arithmetic provider; original v3 fallback and Vulkan property implementation remain, historical property expectations are not a fresh query. C8 already binds actual optional targets.\'\\n"\nsource=source.replace(anchor,extra+anchor,1)\nexec(compile(source,str(basis),\'exec\'))\n'
source=source.replace(anchor,addition,1)
exec(compile(source,'/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/controls/readers/read_controls.py',"exec"))
