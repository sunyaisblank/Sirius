"""Retire exact archived superseded dense providers outside all current closures."""
from pathlib import Path
import hashlib
import sys
R=Path.cwd()
basis=R/'out/retire-native-cpu-predicated-products-provider.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='6a941fea2c59beda683dd9b6109f2a323e84e0b4fef6686e937a947b68b46343'
stage,=sys.argv[1:]
name,bundle={
 'AC':('native-cpu-dense-arithmetic','2e8224da/native-cpu-dense-arithmetic-frame-timeout'),
 'AD':('native-cpu-dense-values','4e49d817/native-cpu-dense-values-frame-timeout')}[stage]
s=basis.read_text()
changes=[
 ("'out/native-cpu-dense-arithmetic-cohort/control-8-acceptance.json'", "'out/native-cpu-masked-arithmetic-cohort/control-8-acceptance.json'"),
 ("('out/native-cpu-dense-arithmetic', 'out/native-cpu-dense-arithmetic-cohort', 'out/native-cpu-dense-arithmetic-compatibility', 'out/native-cpu-dense-values')", "('out/native-cpu-masked-arithmetic', 'out/native-cpu-masked-arithmetic-cohort', 'out/native-cpu-masked-arithmetic-compatibility', 'out/native-cpu-frame-route', 'out/native-cpu-dense-arithmetic', 'out/native-cpu-demand-submit', 'out/native-cpu-demand-submit-cohort')"),
 ('fddb9175/native-cpu-predicated-products-frame-timeout',bundle),
 ('bin/linux-gcc/native-cpu-predicated-products/provider','bin/linux-gcc/'+name+'/provider'),
 ('Two untracked AB provider names only', 'Two untracked '+stage+' provider names only'),
 ('all current AC/contingent AD frozen inputs/build commands/process maps', 'all current AE/route/retained CPP frozen inputs/build commands/process maps'),
 ('AB source/native objects/component/evidence and all active AC/contingent AD inputs/caches/resources retained.', stage+' source/native objects/component/evidence and all active AE/route/retained CPP inputs/caches/resources retained.')]
for a,b in changes:
 assert a in s,a
 s=s.replace(a,b)
exec(compile(s,str(basis),'exec'))
