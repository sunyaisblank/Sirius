"""Retire exact archived superseded raw cohorts, preserving current consumers."""
from pathlib import Path
import hashlib
import sys

R=Path.cwd()
basis=R/'out/retire-native-cpu-predicated-cohort.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest()=='54150716caf5a3490774b75670874c450d2ed013554df907b3a1c7ff5ec6c862'
stage,=sys.argv[1:]
choices={
 'Y':('native-cpu-task-cardinality','bee63496/native-cpu-task-cardinality-frame-timeout'),
 'AA':('native-cpu-static-registers','6c79517c/native-cpu-static-registers-frame-timeout'),
 'AB':('native-cpu-predicated-products','fddb9175/native-cpu-predicated-products-frame-timeout'),
 'AC':('native-cpu-dense-arithmetic','2e8224da/native-cpu-dense-arithmetic-frame-timeout')}
name,bundle=choices[stage]
source=basis.read_text()
changes=[
 ('out/native-cpu-predicated-arithmetic-cohort','out/'+name+'-cohort'),
 ('out/native-cpu-task-cardinality-cohort/control-8-acceptance.json','out/native-cpu-dense-values-cohort/control-8-acceptance.json'),
 ("('out/native-cpu-task-cardinality', 'out/native-cpu-task-cardinality-cohort')",
  "('out/native-cpu-dense-values', 'out/native-cpu-dense-values-cohort', 'out/native-cpu-dense-values-compatibility', 'out/native-cpu-frame-route', 'out/native-cpu-dense-arithmetic')"),
 ('d324ffb6/native-cpu-predicated-arithmetic-frame-timeout',bundle),
 ('raw X trial output outside every current Y input closure.', 'raw '+stage+' trial output outside all current AD/route/retained CPP input closures.'),
 ('raw X cohort files outside every current Y frozen closure','raw '+stage+' cohort files outside every current AD/route/retained CPP frozen closure'),
 ('X native objects and reused full38 standalone component, all current Y inputs/caches/resources preserved.',
  stage+' source/native objects/components and all current AD/route/retained CPP inputs/caches/resources preserved.')]
for before,after in changes:
 assert before in source,before
 source=source.replace(before,after)
exec(compile(source,str(basis),'exec'))
