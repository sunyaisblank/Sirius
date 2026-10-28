from pathlib import Path
import sys
sys.path.insert(0, str(Path.cwd() / 'out/native-cpu-inline-private-cohort'))
from rebind import ao
exec(compile(ao('native_observer.py'), '/home/astra/.project/Sirius/out/native-cpu-v4-cohort/native_observer.py', "exec"))
