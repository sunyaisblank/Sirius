"""Prepare the sole original ABBA only after current C8 and all seven controls."""
from authority import load, an, authority_files, prepare_source
from correctness import check

correctness_inputs = check()
load('prepare.py')
source = an('prepare.py')
assert source.count('out/native-cpu-local-identities/provider') == 2
source = source.replace('out/native-cpu-local-identities/provider', 'out/native-cpu-endpoint-header-provider/provider')
assert source.count('AN local integer identity-store candidate with unchanged AM math') == 1
source = source.replace('AN local integer identity-store candidate with unchanged AM math',
    'AP Endpoint pure/scalar admission extension for header word11241 with six unchanged AO objects including guarded v4/original v3 math')
anchor = 'exec(compile(source, str(basis), "exec"))'
assert source.count(anchor) == 1
source = source.replace(anchor, 'source = prepare_source(source)\n' + anchor, 1)
exec(compile(source, 'out/native-cpu-local-identities-cohort/benefit/prepare.py', 'exec'))
