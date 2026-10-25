"""Reuse validated exact opcode, GOT and MOV-to-LEA target association."""
from pathlib import Path
import hashlib
def associate_fixture(executable,objects):
 R=Path.cwd();base=R/'out/native-cpu-reachable-tag-writers/admit_fixture_build.py'
 source=base.read_text();a=source.index("sys.path.insert(0,str(R/'out/native-cpu-field-compatibility'))");b=source.index("(F/'linked-code-admission.json').write_text",a)
 namespace=dict(R=R,p=dict(executable=str(executable.relative_to(R))),objects=objects)
 exec(compile('from pathlib import Path\nimport re,struct,sys\n'+source[a:b],str(base),'exec'),namespace)
 proof=namespace['proof'];assert proof['passed']
 proof['association_reader']=dict(path=str(base.relative_to(R)),sha256=hashlib.sha256(base.read_bytes()).hexdigest())
 return proof
