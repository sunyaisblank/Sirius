"""Retain assembly labels; complete code equality controls their use."""
from pathlib import Path
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
folder = Path.cwd() / 'bin/linux-gcc/portable-execution-attribution'
output = folder / 'endpoint-labels.o'
assert not output.exists()
subprocess.run(['/usr/bin/as', '--64', '-L', str(folder / 'endpoint.s'), '-o', str(output)], check=True)
assert output.is_file() and output.stat().st_size < 64 * 1024 * 1024
