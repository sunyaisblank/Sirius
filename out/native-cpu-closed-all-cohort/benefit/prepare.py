"""Prepare the sole original ABBA after fresh AX C8 and static composition."""
from authority import load, an, authority_files, prepare_source
from correctness import check

correctness_inputs = check()
load('prepare.py')
source = an('prepare.py')
assert source.count('out/native-cpu-local-identities/provider') == 2
source = source.replace('out/native-cpu-local-identities/provider', 'out/native-cpu-closed-all-provider/provider')
assert source.count('AN local integer identity-store candidate with unchanged AM math') == 1
source = source.replace('AN local integer identity-store candidate with unchanged AM math',
    'AX genuine-job-only closed canonical T/E/DP with existing strict local-effects register math9; unchanged guarded public/scalar/pure admission, ASraw3/GCC/v3/AVX2 fallbacks. Current C8 proves all32 and actual E/DP SIMD/core coverage; Transport12 remains original scalar with no inferred Transport SIMD observation')
anchor = 'exec(compile(source, str(basis), "exec"))'
assert source.count(anchor) == 1
source = source.replace(anchor, 'source = prepare_source(source)\n' + anchor, 1)
exec(compile(source, 'out/native-cpu-local-identities-cohort/benefit/prepare.py', 'exec'))
