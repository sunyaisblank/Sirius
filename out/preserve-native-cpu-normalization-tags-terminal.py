"""Use the frozen preserver with the fourth, fixture-only ISPC source counted."""
from pathlib import Path
import hashlib
R=Path.cwd();p=R/'out/preserve-native-cpu-normalization-tags.py';source=p.read_text()
a="assert len([name for name in source_files if name.endswith('.ispc')]) == 3";b="assert len([name for name in source_files if name.endswith('.ispc')]) == 4"
assert source.count(a)==1;source=source.replace(a,b,1)
a="for path in (Path(__file__).resolve(),):";b="for path in (Path(__file__).resolve(), R/'out/preserve-native-cpu-normalization-tags.py'):"
assert source.count(a)==1;source=source.replace(a,b,1)
exec(compile(source,str(p),'exec'))
