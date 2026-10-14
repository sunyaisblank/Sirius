"""Root-operated provider action; performs no owner-dependent postprocessing."""
from pathlib import Path
import hashlib
import json
import os
import resource
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024**2, 64 * 1024**2))
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
build = root / 'bin/linux-gcc/computed-carry-exit-prerequisites'
producer, mode, folder = sys.argv[1:]
folder = Path(folder).resolve()
assert producer in ('stock', 'candidate') and mode in ('integer', 'normal-transport')
assert folder.parent == work and folder.name == 'provider-' + producer + '-' + mode
fixture = json.loads((work / 'fixture.json').read_text())
for name in ('input', 'expected'):
    assert hashlib.sha256((work / (name + '.bin')).read_bytes()).hexdigest() == fixture[name + '_sha256']
module = build / (producer + '-' + mode + '.spv')
compiled = json.loads((work / 'compile-result.json').read_text())
assert hashlib.sha256(module.read_bytes()).hexdigest() == compiled['modules'][producer + '-' + mode]['sha256']
executable = build / 'gpu'
assert hashlib.sha256(executable.read_bytes()).hexdigest() == compiled['gpu']['sha256']
os.execv(str(executable), [str(executable), str(module), str(work / 'input.bin'),
    str(work / 'expected.bin'), str(folder / 'actual.bin'), mode, '--dispatch'])
