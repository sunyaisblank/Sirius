"""Reuse the exact R full38-case comparator against the newly linked S executable."""
from pathlib import Path
import hashlib
import sys

assert sys.argv[1] == 'run'
old = Path.cwd() / 'out/native-cpu-accumulate-compatibility'
sys.path.insert(0, str(old))
basis = old / 'task.py'
source = basis.read_text()
source = source.replace("'compile-result.json'", "'link-result.json'")
source = source.replace('exact OLD/new raw Transport code', 'exact OLD/new S raw Transport code')
exec(compile(source, str(basis), 'exec'))
