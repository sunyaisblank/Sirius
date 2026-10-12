"""Use the existing LLVM assembler for exact padding/unwind correspondence."""
from pathlib import Path
import json
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
settings = json.loads((work / 'emission-settings.json').read_text())
folder = root / 'bin/linux-gcc/portable-execution-attribution'
output = folder / 'endpoint-llvm-labels.o'
assert not output.exists()
subprocess.run(['/usr/bin/llvm-mc-18', '--filetype=obj', '--save-temp-labels',
                '--large-code-model', '--triple=' + settings['triple'],
                '--mcpu=' + settings['cpu'], '--mattr=' + settings['features'],
                str(folder / 'endpoint.s'), '-o', str(output)], check=True)
assert output.is_file() and output.stat().st_size < 64 * 1024 * 1024
