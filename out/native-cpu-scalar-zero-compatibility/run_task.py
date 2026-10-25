"""Run the unchanged full38 comparator, adapting only receipt/scope labels."""
from pathlib import Path
import sys

assert sys.argv[1] == 'run'
old = Path.cwd() / 'out/native-cpu-accumulate-compatibility'
sys.path.insert(0, str(old))
basis = old / 'task.py'
source = basis.read_text()
source = source.replace("'compile-result.json'", "'build-result.json'")
source = source.replace('all8_native_objects', 'transport_component_objects')
source = source.replace('exact OLD/new raw Transport code', 'exact OLD/new scalar-zero Transport code')
exec(compile(source, str(basis), 'exec'))
