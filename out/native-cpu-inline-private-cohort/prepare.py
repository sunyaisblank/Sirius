"""Keep AO's C8 producer and seal its complete predecessor source chain."""
from rebind import ao, load, authority_files, seal_authorities

ao('prepare.py')
source = load('prepare.py')
source = source.replace('unchanged AM CPP math',
    'same original AM arithmetic with separately guarded optional v4 and original v3 fallback')
anchor = 'exec(compile(source, str(basis), "exec"))'
assert source.count(anchor) == 1
source = source.replace(anchor, 'source = seal_authorities(source)\n' + anchor, 1)
exec(compile(source, '/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/prepare.py', 'exec'))
