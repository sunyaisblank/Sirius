"""Reuse the owned compile loop with exactly the prepared one-object list."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-slang-bridge/build_batch.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == 'b1fec31561d7ae127be4a0ec46a4d036141d4cd80cf14409c36a9b18c921237f'
source = basis.read_text()
source = source.replace('Compile the three sealed typed adapters once; never execute native shader code.',
                        'Compile one pure-admission object once; never execute native shader code.')
source = source.replace("prepared = json.loads((W / 'prepared-batch-inputs.json').read_text())",
    "prepared = json.loads((W / 'prepared-batch-inputs.json').read_text())\n"
    "assert len(prepared['objects']) == 1\nassert not (W / 'batch-compile-result.json').exists()")
source = source.replace('Strict native C++26 compile-only typedABI/packet adapters.',
                        'Strict C++26 compile-only pure-admission C ABI; exactly one new object, three D objects unchanged.')
exec(compile(source, str(basis), 'exec'))
