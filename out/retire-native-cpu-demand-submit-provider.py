"""Retire the exact archived AF provider pair outside all current AG inputs."""
from pathlib import Path
import hashlib
R=Path.cwd();basis=R/'out/retire-native-cpu-predicated-products-provider.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='6a941fea2c59beda683dd9b6109f2a323e84e0b4fef6686e937a947b68b46343'
s=basis.read_text()
changes=[
 ("'out/native-cpu-dense-arithmetic-cohort/control-8-acceptance.json'", "'out/native-cpu-arithmetic-closure-cohort/control-8-acceptance.json'"),
 ("('out/native-cpu-dense-arithmetic', 'out/native-cpu-dense-arithmetic-cohort', 'out/native-cpu-dense-arithmetic-compatibility', 'out/native-cpu-dense-values')", "('out/native-cpu-arithmetic-closure', 'out/native-cpu-arithmetic-closure-cohort', 'out/native-cpu-arithmetic-closure-compatibility', 'out/native-cpu-masked-arithmetic', 'out/native-cpu-frame-route', 'out/native-cpu-dense-arithmetic')"),
 ('fddb9175/native-cpu-predicated-products-frame-timeout','6176354b/native-cpu-demand-submit-frame-timeout'),
 ('bin/linux-gcc/native-cpu-predicated-products/provider','bin/linux-gcc/native-cpu-demand-submit/provider'),
 ('Two untracked AB provider names only','Two untracked superseded AF provider names only'),
 ('all current AC/contingent AD frozen inputs/build commands/process maps','all current AG/retained arithmetic/route frozen inputs/build commands/process maps'),
 ('AB source/native objects/component/evidence and all active AC/contingent AD inputs/caches/resources retained.','AF source/native objects/evidence and all active AG/retained arithmetic/route inputs/caches/resources retained.')]
for a,b in changes:
 assert a in s,a
 s=s.replace(a,b)
exec(compile(s,str(basis),'exec'))
