"""Prepare the one raw comparison after completed read-only build accounting."""
from pathlib import Path
import sys

assert sys.argv[1:] == []
sys.argv.append('run')
basis = Path(__file__).resolve().parent / 'prepare.py'
source = basis.read_text()
assert source.count("owner_ok(W / 'compile/owner.json')") == 1
source = source.replace("owner_ok(W / 'compile/owner.json')", "owner_ok(W / 'proof/owner.json')")
source = source.replace("'compile-result.json'", "'compile-completion.json'")
source = source.replace("'compile/owner.json'", "'proof/owner.json'")
source = source.replace("str(W / 'task.py')", "str(W / 'run_task.py')")
source = source.replace("inputs.update(build['runtime_libraries'])", "inputs.update(build['runtime_libraries'])\n    inputs.update(build['proof_completion']['source_inputs'])")
exec(compile(source, str(basis), 'exec'))
