"""Prepare the sole original ABBA after fresh AR C8 and static composition."""
from authority import load, an, authority_files, prepare_source
from correctness import check

correctness_inputs = check()
load('prepare.py')
source = an('prepare.py')
assert source.count('out/native-cpu-local-identities/provider') == 2
source = source.replace('out/native-cpu-local-identities/provider', 'out/native-cpu-packed-math-provider/provider')
assert source.count('AN local integer identity-store candidate with unchanged AM math') == 1
source = source.replace('AN local integer identity-store candidate with unchanged AM math',
    'AR private math-boundary raw3 and reviewed zero-residual optional arithmetic; unchanged AP first7 objects, baseline selector and original v3 fallback')
anchor = 'exec(compile(source, str(basis), "exec"))'
assert source.count(anchor) == 1
source = source.replace(anchor, 'source = prepare_source(source)\n' + anchor, 1)
exec(compile(source, 'out/native-cpu-local-identities-cohort/benefit/prepare.py', 'exec'))
