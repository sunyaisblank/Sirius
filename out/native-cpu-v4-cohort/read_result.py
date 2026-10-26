"""Keep every original C8 gate and add scoped optional-target joins."""
from pathlib import Path
from rebind import load

source = load('read_result.py')
anchor = "exec(compile(source, str(basis), 'exec'))"
assert source.count(anchor) == 1
addition = '''
anchor = "(work / 'control-8-acceptance.json').write_text"
assert source.count(anchor) == 1
extra = """
from optional_math import check as check_optional_math
result['optional_math'] = check_optional_math(root, action, observation)
result['hardware_expectations_scope'] = 'Historical providerA float controls/features and same-CPU identity supply expectations only. AO changes only private object composition/identity, preserving actual existing v3 eligibility, FP/precision properties, raw3/scheduler/fallback and original physical contracts. New optional v4 uses exact original arithmetic and separately guarded ISA. Current C8 observes selected mode, native rows/pure/SIMD/scalar returns and bounded optional PCs/cachebit. No fresh property query or older-hardware execution is claimed.'
result['scope'] = 'Sole original C8; original51markers54modulejoins32complete-bufferdigests, actual T/E/DP rows/pure returns and selected original scheduling paths. Optional math runtime coverage is separately explicit. No benefit/all8/frame/scientific/release/adoption claim.'
"""
source = source.replace(anchor, extra + '\\n' + anchor, 1)
exec(compile(source, str(basis), 'exec'))
'''
source = source.replace(anchor, addition, 1)
exec(compile(source, str(Path.cwd() / 'out/native-cpu-local-identities-cohort/read_result.py'), 'exec'))
