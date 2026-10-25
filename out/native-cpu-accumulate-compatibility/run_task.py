"""Execute the one original comparison using the completed single-compile build."""
from pathlib import Path
import sys

assert sys.argv[1] == 'run'
basis = Path(__file__).resolve().parent / 'task.py'
source = basis.read_text().replace("'compile-result.json'", "'compile-completion.json'")
exec(compile(source, str(basis), 'exec'))
