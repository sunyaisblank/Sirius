from rebind import load
source=load('prepare.py')
source=source.replace('unchanged AM CPP math','same original AM arithmetic with separately guarded optional v4 and original v3 fallback')
exec(compile(source, '/home/astra/.project/Sirius/out/native-cpu-local-identities-cohort/prepare.py', "exec"))
