"""Original sole cold frame, H products/resources, AO only after all8/all22."""
from pathlib import Path
from frame_authority import load
source=load('prepare.py')
anchor="assert c8['passed'] and c8['all32_readbacks_match_preserved_originals']"
assert source.count(anchor)==1
source=source.replace(anchor,anchor+"\nassert c8['optional_math']['all_selected_optional_math_observed']",1)
before='exact AN provider with original AM identity bank stores before C++ calls, unchanged general arithmetic/AK refinements and scheduler/FP scope'
after='exact AO provider with separately guarded optional v4 arithmetic and byte-preserved original v3 fallback, original AM identities before C++ calls and unchanged original general expressions/AK refinements/raw3/scalar tails/scheduler/FP scope'
assert source.count(before)==1
source=source.replace(before,after,1)
exec(compile(source,'/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/cold-frame/prepare.py',"exec"))
